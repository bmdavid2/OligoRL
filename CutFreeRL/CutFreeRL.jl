

import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences, XLSX,ArgParse

# Overload the concatenation operator to combine a sequence and 
# a single base, i.e. dna"AGCGTGC" * DNA_T
function *(seq::LongSequence{DNAAlphabet{4}}, base::DNA)
    seq * LongDNA{4}([base])
end

# Overload sort to order DNA codes by degeneracy:
#   A,G,C,T,M,R,W,S,Y,K,V,H,D,B,N
function sort(seq::LongSequence{DNAAlphabet{4}}; rev::Bool=false)
    sort(collect(seq), by=degeneracy, rev=rev)
end

"""
    isvalid(seq, sites)

Check if any of the sites appear in the sequence.
"""
function isvalid(seq, sites)
    for site in sites
        if occursin(ExactSearchQuery(site,iscompatible), seq)
            return false
        end
    end
    return true
end

"""
    find_valid_randomer(prefix, bases, sites)

Randomly select bases to add the prefix, ensuring no site
appears in the randomer.
"""
function find_valid_randomer(prefix, bases, sites)
    n = length(bases)
    pre_n = length(prefix)
    randomer = prefix * (dna"-" ^ n)
    for i = 1:n
        candidates = shuffle(bases[i])
        for cand in candidates
            randomer[pre_n+i] = cand
            if isvalid(randomer[1:pre_n+i], sites)
                break
            else
                randomer[pre_n+i] = DNA_Gap
            end
        end
        if randomer[pre_n+i] == DNA_Gap
            break
        end
    end
    return randomer
end

"""
    degeneracy(base::DNA[, uselog2=true])

Calculate the degeneracy of a single base. The degeneracy 
is the number of bases in the code, so
    degeneracy(DNA_N) == 4
and
    degeneracy(DNA_A) == 1

If `uselog2`, the degeneracy is log2 transformed. Because degeneracy(DNA_Gap)=0, we approximate
the log2 transformation as 2^-100. 
"""
function degeneracy(base::DNA; uselog2=true)
    deg = 0.0
    if base == DNA_N
        deg = 4.0
    elseif base == DNA_A || base == DNA_C || base == DNA_G || base == DNA_T
        deg = 1.0
    elseif base == DNA_B || base == DNA_D || base == DNA_H || base == DNA_V
        deg = 3.0
    elseif (base == DNA_K || base == DNA_Y || base == DNA_S || base == DNA_W 
            || base == DNA_R || base == DNA_M)
        deg = 2.0
    else # DNA_Gap
        deg = 0.0
    end
    
    if uselog2
        if base == DNA_Gap
            deg = 2^-100
        end
        return log2(deg)
    else
        return deg
    end
end

"""
    degeneracy(seq::LongSequence{DNAAlphabet{4}}; uselog2=true)

Calculate the degeneracy of a DNA sequence.
"""
function degeneracy(seq::LongSequence{DNAAlphabet{4}}; uselog2=true)
    if uselog2
        return sum(degeneracy.(seq, uselog2=true))
    else
        return prod(degeneracy.(seq ;uselog2=false))
    end
end

# Run Monte Carlo simulations using a random base policy over a fixed horizon.
function simulate_random(prefix, bases, sites; nsims=1000, horizon=length(bases))
    degs = zeros(nsims)
    true_horizon = min(horizon, length(bases)) # don't go beyond the randomer length
    for i = 1:nsims
        degs[i] = degeneracy(find_valid_randomer(prefix, bases[1:true_horizon], sites))
    end
    return mean(degs)
end

# Run a 1-step greedy lookahead policy.
function simulate_greedy(prefix, bases, sites; horizon=length(bases))
    true_horizon = min(horizon, length(bases)) # don't go beyond the randomer length
    for h = 1:true_horizon
        # Iterate through the bases by decending degeneracy, so the greedy approach
        # is to stop as soon as we've found one.
        found = false # define out here for scoping rules
        for base in sort(bases[h], rev=true)
            found = false
            if isvalid(prefix * base, sites)
                prefix *= base
                found = true
                break
            end
        end
        if ~found
            # None of the bases were valid, so stop and leave the DNA_Gap
            break
        end
    end
    return degeneracy(prefix)
end

"""
    CutFreeRL(sequence::LongSequence{DNAAlphabet{4}}, sites::Array{LongSequence{DNAAlphabet{4}},1}; simulate=simulate_random, kwargs...)

Run the rollout algorithm to solve the CutFree MDP. 

# Arguments 
- `sequence`: Starting DNA sequence that should be blocked from containing restriction sites. To generate a set of barcodes with the highest diversity, start with a string of N's the length of your oligo. Option is given because some companies restrict which degernate bases are allowed.
- `sites`: An array of restriciton enzyme recognition sequences to be blocked in the random barcode. 
- `simulate`: The choice of policy for rollout simulations. simulate_random will use a random rollout policy. simulate_greedy will use a greedy 1 step lookahead policy. 
#Optional Keyword Arugments 
- `nsims`: The number of rollout simulations per action
"""
function CutFreeRL(sequence::LongSequence{DNAAlphabet{4}}, sites::Array{LongSequence{DNAAlphabet{4}},1}; simulate=simulate_random, kwargs...)
    allowedbasedict=Dict(DNA_A => dna"A", DNA_T => dna"T", DNA_C => dna"C", DNA_G => dna"G",DNA_R =>dna"AG",DNA_Y=>dna"CT",DNA_S=>dna"GC",
    DNA_W=>dna"AT",DNA_K=>dna"GT",DNA_M=>dna"AC",DNA_B=>dna"CGTYSKB",DNA_D=>dna"AGTRWKD",DNA_H=>dna"ACTYWMH",DNA_V=>dna"ACGRSMV",DNA_N=>dna"ACGTRYSWKMBDHVN",DNA_Gap=>dna"-")
    bases=[dna"-" for i=1:length(sequence)]
    for i in eachindex(sequence)
        bases[i]=allowedbasedict[sequence[i]]
    end 
    n = length(bases)
    randomer = dna"-" ^ n

    # If the restriction site is non-palindromic, we also need to block
    # the reverse complement.
    for i = 1:length(sites)
        if ~ispalindromic(sites[i])
            push!(sites, reverse_complement(sites[i]))
        end
    end

    # Find the longest restriction site and use this as the horizon.
    max_len = map(length, sites) |> maximum
    horizon = max_len

    for i = 1:n
        best_base = DNA_Gap
        best_deg = -1.0
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        for base in bases[i]
            if ~isvalid(randomer[start:i-1] * base, sites)
                # adding this bases creates a restriction site; skip it
                continue
            end
            deg = simulate(randomer[start:i-1] * base, bases[i+1:end], sites; horizon=horizon, kwargs...)

            # Update use this base if:
            #   1. the mean randomer degeneracy is higher than the previous best, OR
            #   2. the randomer degeneracy is tied with the previous best, but the new base has
            #      higher individual degeneracy
            if deg > best_deg || (deg == best_deg && degeneracy(base) > degeneracy(best_base))
                best_deg = deg
                best_base = base
            end
        end

        if best_base == DNA_Gap
            # The sequence has terminated; there are no bases that can be added without making
            # a restriction site.
            break
        else
            randomer[i] = best_base
        end
    end

    return randomer
end

#= Test code
sequence=dna"NNNNNNNNNNNNNNNNNNNN"
sites = [dna"GGTCTC"] # BsaI
randomer = CutFreeRL(sequence, sites, simulate=simulate_random, nsims=1000)
=#

#
#randomer = CutFreeRL(sequence, sites, simulate=simulate_greedy)

#deg = degeneracy(randomer)
#natdeg = deg / log2(exp(1))
#println("Final randomer: $randomer")
#println("Final degeneracy: $deg ($natdeg)")

# Run the benchmark CutFree set
function runtimes(data)
    data[:random_objval] = 0.0
    data[:random_runtime] = 0.0
    data[:random_code] = ""

    data[:greedy_objval] = 0.0
    data[:greedy_runtime] = 0.0
    data[:greedy_code] = ""

    #data.objval = data.objval ./ log(2)

    #bases = [all_bases for _ in 1:20] #Used for running sites data 
    sites = map(x -> split(x, ',') .|> LongDNASeq, data.sites)
    for i = 1:nrow(data)
        bases = [all_bases for _ in 1:data.oligo_lengths[i]] #used for running lengths data 
        stats = @timed CutFreeRL(bases, sites[i], nsims=1000)
        data[i,:random_objval] = degeneracy(stats.value; uselog2=false)
        data[i,:random_runtime] = stats.time
        data[i,:random_code] = convert(String, stats.value)

        stats = @timed CutFreeRL(bases, sites[i], simulate=simulate_greedy)
        data[i,:greedy_objval] = degeneracy(stats.value;uselog2=false)
        data[i,:greedy_runtime] = stats.time
        data[i,:greedy_code] = convert(String, stats.value)
        println(i)
    end
    return data
end

#data = runtimes(CSV.read("./results.csv")[1:10,:])
#data = runtimes(CSV.read("./results_sites_4_1_21.csv"))
#CSV.write("rollout_results_sites_4_1_21.csv",data)
#data = runtimes(CSV.read("./results_length_4_15_21_combined.csv"))
#CSV.write("rollout_results_length_4_15_21.csv",data)
function read_restriction_sites(file="cutfree_rebase_data.xlsx") #File taken from the rebase.R code. Filtered for the 56 commerially available, palindromic, 6 bp, restriction sites 
    xf=XLSX.readxlsx(file);
    x=xf[XLSX.sheetnames(xf)[1]]
    blocking_sites=[]
    dims=size(x[:]);
    n=dims[1]
 for i in 2:n #start from second row because excel uses row 1 as the label 
        seq=LongDNA{4}(string(x[i,7])) #the sequences are the 7th column of this dataframe
        push!(blocking_sites,seq)
    return blocking_sites
end
end


function random_subset(data;size=10)
    shuffled=shuffle(data)
    if size >length(data)
        return shuffled
    else 
        return shuffled[1:size]
    end 
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--sequence","-s"
            help="Starting DNA sequence that should be blocked from containing restriction sites. To generate a set of barcodes with the highest diversity, start with a string of N's the length of your oligo. "
            required=true   
        "--restrictionsites","-r"
            help = "Sequences to block from the oligo pools. Separate multiple sequences by commas. Do not include spaces."
            required=true
        "--nsims","-n"
            help = "Number of simulations per rollout"
            arg_type= Int
            default=100
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end


    sequence=parsed_args["sequence"]
    sites=parsed_args["restrictionsites"]
    nsims=parsed_args["nsims"]
    sequence=LongDNA{4}.(sequence)
    sites=LongDNA{4}.(split(sites,","))
    randomer=CutFreeRL(sequence,sites;nsims=nsims)
    randomer=String(randomer)
    show(stdout,"text/plain",randomer)
    println("\n")
end


main()


