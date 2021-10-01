import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences, FASTX, XLSX, CPUTime
import BioSequences: iscompatible



function *(seq::LongDNASeq, base::DNA)
    seq * LongDNASeq([base])
end

# Overload sort to order DNA codes by degeneracy:
#   A,G,C,T,M,R,W,S,Y,K,V,H,D,B,N
function sort(seq::LongSequence{DNAAlphabet{4}}; rev::Bool=false)
    sort(convert(Vector{DNA}, seq), by=degeneracy, rev=rev)
end

# Determine whether or not a candidate sequence belongs to any member of the list of sites. We only want seq that belong to sites 
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
            deg = 2^-100  #log2(0) is negative infinity, so we use something close to zero 
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
        return prod(degeneracy.(seq, uselog2=false))
    end
end

"""
    pooled_degeneracy(pool)

Calcualte the sum degeneracy of pool.
"""
function pooled_degeneracy(pool)
    deg=0;
    for randomer in pool
        deg +=degeneracy(randomer;uselog2=false);
    end
    return deg
end 




"""
    iscompatible(seq_A::LongSequence{DNAAlphabet{4}},seq_B::LongSequence{DNAAlphabet{4}})

Check if seq_A and seq_B contain at least one common sequence. Seq_A and Seq_B must be the same length.
"""
function iscompatible(seq_A::LongSequence{DNAAlphabet{4}},seq_B::LongSequence{DNAAlphabet{4}})
    if length(seq_A) == length(seq_B)
        return prod(iscompatible.(seq_A,seq_B))
    else # Sequences are different lenghts 
        return false 
    end 
end 

###############
#OligoRL algorithm
###############



# Run Monte Carlo simulations using a random base policy over a fixed horizon.
function simulate_random(prefix, bases, sites; nsims=1000, horizon=length(bases), reward=degeneracy)
    sim_reward = zeros(nsims)
    true_horizon = min(horizon, length(bases)) # don't go beyond the randomer length
    for i = 1:nsims
        sim_reward[i] = reward(find_valid_randomer(prefix, bases[1:true_horizon], sites))
    end
    return mean(sim_reward)
end

# Run a 1-step greedy lookahead policy.
function simulate_greedy(prefix, bases, sites; horizon=length(bases), reward=degeneracy)
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
    return reward(prefix)
end

function rollout(bases::Array{LongSequence{DNAAlphabet{4}},1}, sites::Array{LongSequence{DNAAlphabet{4}},1}; simulate=simulate_random, kwargs...)
    n = length(bases)
    randomer = dna"-" ^ n

    # If the unwanted site is non-palindromic, we also need to block
    # the reverse complement.
    for i = 1:length(sites)
        if ~ispalindromic(sites[i])
            push!(sites, reverse_complement(sites[i]))
        end
    end

    # Find the longest unwanted site and use this as the horizon.
    max_len = map(length, sites) |> maximum
    horizon = max_len

    for i = 1:n
        best_base = DNA_Gap
        best_score = -1.0
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        for base in bases[i]
            if ~isvalid(randomer[start:i-1] * base, sites)
                # adding this bases creates a restriction site; skip it
                continue
            end
            score = simulate(randomer[start:i-1] * base, bases[i+1:end], sites; horizon=horizon, kwargs...)

            # Update use this base if:
            #   1. the mean randomer score is higher than the previous best
            if score > best_score 
                best_score = score
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


randomerlen=10
all_ns=dna"AGCTMRWSYKVHDBN"
all_bases=[ all_ns for i =1:randomerlen]

