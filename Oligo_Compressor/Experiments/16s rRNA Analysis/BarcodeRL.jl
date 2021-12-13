

import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences, Combinatorics

# Overload the concatenation operator to combine a sequence and 
# a single base, i.e. dna"AGCGTGC" * DNA_T
function *(seq::LongDNASeq, base::DNA)
    seq * LongDNASeq([base])
end

# Overload sort to order DNA codes by degeneracy:
#   A,G,C,T,M,R,W,S,Y,K,V,H,D,B,N
function sort(seq::LongSequence{DNAAlphabet{4}}; rev::Bool=false)
    sort(convert(Vector{DNA}, seq), by=degeneracy, rev=rev)
end


function is_palindrome(randomer)
    return occursin(randomer,reverse_complement(randomer))
end


"""
    isvalid(seq, sites)

Check if any of the sites appear in the sequence.
"""
function isvalid(seq, sites)
    for site in sites
        if occursin(site, seq)
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

If `uselog2`, the degeneracy is log2 transformed. By convention,
the degneracy of `DNA_Gap` is zero in both the log2-transformed
and the untransformed case.
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


## Read in the list of sequences you want to exclude in the designed primers. 
# Function finds every unique hexamer in the list of sequences and stores it as an array
#. For example, use list of rRNA and tRNA sequences. This funciton takes care of orientation. We end up using the reverse complement since we are finding primers that won't hybridize 
## Example
#blocking_sites=read_blocking_sites(CSV.read("./SMU_UA159_rRNA_tRNA.csv"))
function read_invalid_sites(data; randomerlen=6)
    blocking_sites=[];
    for i =1:nrow(data)
        seq=reverse_complement(LongDNASeq(string(data[i,:Sequence])));
        for j=1:length(seq)-(randomerlen-1)
            if isvalid(seq[j:j+(randomerlen-1)],blocking_sites)
                push!(blocking_sites,seq[j:j+(randomerlen-1)])
            end
        end
    end 
    return blocking_sites
end 


function edit_distance(seq1::LongSequence{DNAAlphabet{4}},seq2::LongSequence{DNAAlphabet{4}})
    m=length(seq1);
    n=length(seq2);
    d=zeros(m+1,n+1)
    for i = 1:m
        d[i+1,1]=i
    end 
    for j = 1:n
        d[1,j+1]=j
    end 
    for j = 1:n
        for i=1:m
            substitution_cost=1
            if iscompatible(seq1[i],seq2[j])
                substitution_cost=0;
            end 
            d[i+1,j+1]= min(d[i,j+1]+1,d[i+1,j]+1,d[i,j]+substitution_cost)
        end 
    end 
    return d[m+1,n+1]
end 

function hamming_distance(seq1::LongSequence{DNAAlphabet{4}},seq2::LongSequence{DNAAlphabet{4}})
    dist_counter=0
    if  length(seq1) != length(seq2)
        error("Sequences must be the same length.")
    end 
    n=length(seq1)
    for i=1:n
        if !iscompatible(seq1[i],seq2[i])
            dist_counter+=1
        end 
    end 
    return dist_counter
end 
function reward_function_factory(barcode_pool)
    function reward_function(barcode;distance=hamming_distance)
        sum_distances=0
        min_distance=length(barcode)
        for code in barcode_pool
            dist=distance(barcode,code)
            sum_distances+=dist
            min_distance=min(min_distance,dist)

        end
        if is_palindrome(barcode)
            sum_distances=0;
            min_distance=0;
        end 
        return mean(sum_distances), min_distance
    end 
    return reward_function
end 




function generate_available_ns(base_randomer::LongSequence{DNAAlphabet{4}})
    randomerlen=length(base_randomer)
    basedict=Dict(DNA_A => dna"A", DNA_T => dna"T", DNA_C => dna"C", DNA_G => dna"G",DNA_R =>dna"AG",DNA_Y=>dna"CT",DNA_S=>dna"GC",
    DNA_W=>dna"AT",DNA_K=>dna"GT",DNA_M=>dna"AC",DNA_B=>dna"CGT",DNA_D=>dna"AGT",DNA_H=>dna"ACT",DNA_V=>dna"ACG",DNA_N=>dna"ACGT",DNA_Gap=>dna"-")

    available_ns=[dna"-" for _ in 1:randomerlen]
    for i=1:randomerlen
        available_ns[i]=basedict[base_randomer[i]]
    end 
    return available_ns
end 

function generate_randomer_sample(available_ns)
    randomerlen=length(available_ns);
    randomer=dna"-"^randomerlen
    for i =1:randomerlen
        randomer[i]=shuffle(available_ns[i])[1]
    end 
    return randomer
end 


function simulate_random(prefix, bases, sites; nsims=1000, horizon=length(bases),reward=reward_function,kwargs...)
    score = zeros(nsims)
    min_distance=zeros(nsims)
    true_horizon = min(horizon, length(bases)) # don't go beyond the randomer length
    for i = 1:nsims
        score[i],min_distance[i] = reward(find_valid_randomer(prefix, bases[1:true_horizon], sites))
    end
    return mean(score),mean(min_distance)
end





function oligo_rollout(bases, sites; simulate=simulate_random, kwargs...)
    n = length(bases)
    randomer = dna"-" ^ n


    # Find the longest restriction site and use this as the horizon.
    max_len = map(length, sites) |> maximum
    horizon = max_len

    for i = 1:n
        best_base = DNA_Gap
        best_min_distance = 0
        best_score=0
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        for base in bases[i]
            if ~isvalid(randomer[start:i-1] * base, sites)
                # adding this bases creates a restriction site; skip it
                continue
            end
            score,min_distance = simulate(randomer[start:i-1] * base, bases[i+1:end], sites; horizon=horizon, kwargs...)

            # Update use this base if:
            #   1. the mean randomer degeneracy is higher than the previous best, OR
            #   2. the randomer degeneracy is tied with the previous best, but the new base has
            #      higher individual degeneracy
            if min_distance > best_min_distance || (min_distance== best_min_distance && score > best_score)
                best_min_distance=min_distance
                best_score= score
                best_base = base
            end
        end

        if best_base == DNA_Gap
            # The sequence has terminated; there are no bases that can be added without making
            # a valid barcode
            break
        else
            randomer[i] = best_base
        end
    end

    return randomer
end


function BarcodeRL(sites;nbarcodes=16,barcodelen=8,kwargs...)
    barcode_pool=[];
    base_barcode_ns=dna"N"^barcodelen
    available_ns=generate_available_ns(base_barcode_ns)
    for i= 1:nbarcodes
        reward_function=reward_function_factory(barcode_pool)
        barcode=oligo_rollout(available_ns,sites;reward=reward_function, kwargs...)
        push!(barcode_pool,barcode)
        push!(sites,barcode)
        push!(sites,reverse_complement(barcode))
    end
    reward_function=reward_function_factory(barcode_pool[1:end-1])
    score,min_distance=reward_function(barcode_pool[end];distance=hamming_distance)
    return barcode_pool,score,min_distance
end 

function generate_initial_population(nbarcodes,barcode_form,population_size)
    barcodelen=length(barcode_form)
    initial_population=[[dna"-"^barcodelen for i =1:nbarcodes] for j=1:population_size]
    available_ns=generate_available_ns(barcode_form)
    for j=1:population_size
        for k=1:nbarcodes 
            for i=1:barcodelen
                initial_population[j][k][i]=shuffle(available_ns[i])[1]
            end 
        end 
    end 
    return initial_population

end 

function fitness(barcodes,sites,x_combs)
    n=length(barcodes)
    fits=zeros(length(x_combs))
    for i=1:length(x_combs)
        fits[i]=hamming_distance(barcodes[x_combs[i][1]],barcodes[x_combs[i][2]])
    end 
    fitness=minimum(fits)
    for i=1:length(barcodes)
        if ~isvalid(barcodes[i],vcat(sites,barcodes[union(1:i-1,i+1:n)]))
            fitness=0
            break
        end 
    end 
    return fitness
end 

function mutation(barcodes,sites,barcode_form;mutation_prob=0.25,kwargs...)
    new_barcodes=deepcopy(barcodes)
    n=length(new_barcodes)
    mutation_sites=findall(y-> y==1,rand(n).<mutation_prob)
    available_ns=generate_available_ns(barcode_form)
    order=shuffle(mutation_sites)
    for i in order
        barcode_pool=new_barcodes[union(1:i-1,i+1:n)]
        reward_function=reward_function_factory(barcode_pool)
        new_barcodes[i]=oligo_rollout(available_ns,sites;reward=reward_function,kwargs...)
    end 
    return new_barcodes
end 
#barcodes=[dna"ATGCGT",dna"ATGGGG",dna"GCACGT",dna"ATTGAT"]
#sites=read_invalid_sites(CSV.read("./Invalid_Sequences_Test.csv",DataFrame);randomerlen=6)
#test=mutation(barcodes,sites)

function pair_off(candidates)
    npairs=Int(floor(length(candidates)/2))
    pairs=[zeros(Int64,2) for i =1:npairs]
    cands=shuffle(candidates)
    for i =1:npairs
        pairs[i][1]=cands[2*i-1]
        pairs[i][2]=cands[2*i]
    end 
    return pairs

end 
function recombination(population;recombination_prob=0.25,kwargs...)
    new_population=deepcopy(population)
    n_barcodes=length(population[1])
    n_chromosomes=length(population);
    recombination_sites=findall(y-> y==1, rand(n_chromosomes).< recombination_prob)
    r_pairs=pair_off(recombination_sites)
    for i=1:length(r_pairs)
        r_index=Int(ceil(rand()*n_barcodes))
        new_population[r_pairs[i][1]][r_index:end],new_population[r_pairs[i][2]][r_index:end]=new_population[r_pairs[i][2]][r_index:end],new_population[r_pairs[i][1]][r_index:end];
    end 

    return new_population[recombination_sites]
end
#Test 
#population=[[dna"AAAA",dna"AAAA",dna"AAAA",dna"AAAA"],[dna"GGGG",dna"GGGG",dna"GGGG",dna"GGGG"]]
#recombination(population;recombination_prob=1)

function rank_probability(n_chromosomes;Pc=0.25,kwargs...)
    probs=zeros(n_chromosomes)
    cumul_prob=zeros(n_chromosomes)
    for i=1:n_chromosomes-1
        probs[i]=Pc*(1-Pc)^(i-1)
        cumul_prob[i]=sum(probs[1:i])
    end 
    probs[end]=(1-Pc)^(n_chromosomes-1)
    cumul_prob[end]=sum(probs)
    return cumul_prob
end 


function selection(population,fitnesses;kwargs...)
    n_chromosomes=length(population)
    ranks=sortperm(fitnesses,rev=true)
    cumul_prob=rank_probability(n_chromosomes;kwargs...)
    sample=rand()
    cutoff=minimum(findall(y-> y>=sample,cumul_prob));
    idx=ranks[cutoff]
    return population[idx]
end
#Test
#population=[[dna"AAAA",dna"AAAA",dna"AAAA",dna"AAAA"],[dna"GGGG",dna"GGGG",dna"GGGG",dna"GGGG"]]
#fitnesses=[20, 10]
#selection(population,fitnesses)

#results=zeros(Int64,100000)
#for i=1:100000
#    results[i]=selection(population,fitnesses)
#end 
#x=length(findall(y->y==1,results))
#t=length(findall(y->y==2,results))
#thy=length(findall(y->y==30,results))
#thy/100000

    





function GA_BarcodeRL(sites;nbarcodes=16,barcode_form=dna"NNNNNNNN",generations=50,population_size=30,kwargs...)
    x_comb=collect(combinations(1:nbarcodes,2));
    population=generate_initial_population(nbarcodes,barcode_form,population_size)
    for i=1:generations
        for j=1:population_size
            population=push!(population,mutation(population[j],sites,barcode_form;kwargs...))
        end
        recomb_pop=recombination(population;kwargs...)
        for j =1:population_size-1
            recomb_pop=vcat(recomb_pop,recombination(population;kwargs...))
        end 
        population=vcat(population,recomb_pop)
        fitnesses=fitness.(population,(sites,),(x_comb,))
        new_pop_indicies=zeros(Int64,population_size)
        v_indicies=collect(1:length(population))
        for i=1:population_size
            new_pop_cand=selection(population[v_indicies],fitnesses[v_indicies];kwargs...)
            new_pop_indicies[i]=findall(y-> y == new_pop_cand,population)[1]
            v_indicies=filter(y->y != new_pop_indicies[i],v_indicies)
        end 
        population=population[new_pop_indicies]
        println("Generation: $i")
    end 

    final_fitness=fitness.(population,(sites,),(x_comb,))
    max_fitness,max_pop_idx=findmax(final_fitness)
    return population[max_pop_idx], max_fitness
end 





sites=read_invalid_sites(CSV.read("./Invalid_Sequences_Test.csv",DataFrame);randomerlen=8)

barcodes,best_fitness=GA_BarcodeRL(sites;nbarcodes=16,barcode_form=dna"NNNNNNNN",generations=10,population_size=30,nsims=100,Pc=0.25,recombination_prob=0.25,mutation_prob=0.25)

function benchmark_GA_barcodeRL()
    num_barcodes=repeat([4 ,8 ,12, 16, 20, 24, 28, 32, 36, 40, 44, 48],3)
    nexps=length(num_barcodes)
    time=zeros(nexps)
    min_edit_distance=zeros(nexps)
    sites=read_invalid_sites(CSV.read("./Invalid_Sequences_Test.csv",DataFrame);randomerlen=8)
    barcodes=["-" for i=1:nexps]
    for i =1:nexps
        data= @timed GA_BarcodeRL(sites;nbarcodes=num_barcodes[i],barcode_form=dna"NNNNNNNN",generations=10,population_size=50,nsims=100,Pc=0.25,recombination_prob=0.25,mutation_prob=0.25)
        time[i]=data.time
        min_edit_distance[i]=data.value[2]
        barcodes[i]=join(convert.(String,data.value[1]),",")
    end 
    DF=DataFrame(num_barcodes=num_barcodes,barcodes=barcodes,min_edit_distance=min_edit_distance,time=time)
    return DF
end 


data=benchmark_GA_barcodeRL()
CSV.write("GA_BarcodeRL_benchmark_replicates.csv",data)









