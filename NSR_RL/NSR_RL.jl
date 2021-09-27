import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences
import BioSequences: iscompatible

################################################
#OligoRL Toolset 
#################################################

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


# Determine if a given sequence is palindromic. Works both for unambiguous and ambigous bases. Returns TRUE if the sequence is a palindrome. 

"""
    is_palindrome(seq)

Check if seq is palindromic. 
"""
function is_palindrome(seq)
    return occursin(seq,reverse_complement(seq))
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

#Calculate the degeneracy of an entire oligo pool

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

#Determine if two oligos contain at least one common sequence. Used for assessing ambigous oligos. 
#Extends the iscompatible function in the BioSequences package from bases to sequences.
#The two sequences must be the same length, otherwise the function returns false 

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

#Decompress an ambiguous oligo into all of its possible sequences. Input is a single oligo.
#
"""
    decompress_oligo(comp_oligo)

Expand a degnerate oligo into the non-degnerate oligos it defines. 
"""

function decompress_oligo(compressed_oligo::LongSequence{DNAAlphabet{4}})
    basedict=Dict(DNA_A => dna"A", DNA_T => dna"T", DNA_C => dna"C", DNA_G => dna"G",DNA_R =>dna"AG",DNA_Y=>dna"CT",DNA_S=>dna"GC",
    DNA_W=>dna"AT",DNA_K=>dna"GT",DNA_M=>dna"AC",DNA_B=>dna"CGT",DNA_D=>dna"AGT",DNA_H=>dna"ACT",DNA_V=>dna"ACG",DNA_N=>dna"ACGT",DNA_Gap=>dna"-")
    nseqs=trunc(Int128,degeneracy(compressed_oligo;uselog2=false))
    randomerlen=length(compressed_oligo)
    blockarray=zeros(Int128,randomerlen+1)
    blockarray[1]=nseqs
    degsarray=zeros(Int128,randomerlen)
    repsarray=zeros(Int128,randomerlen+1)
    repsarray[1]=1
    uncompressed_oligos=[dna"-"^randomerlen for _ in 1:nseqs]
    for i = 1:randomerlen 
        nucleotides=basedict[compressed_oligo[i]]
        degsarray[i]=length(nucleotides)
        blockarray[i+1]=trunc(Int128,blockarray[i]/degsarray[i]);
        bases=dna"-"^blockarray[i];
        repsarray[i+1]=trunc(Int128,prod(degsarray[1:i]));
            for  k in 1:degsarray[i]
                for l in 1:blockarray[i+1]
                    bases[l+blockarray[i+1]*(k-1)]=nucleotides[k]
                end
            end
        
        bases=bases^repsarray[i];
        for m in 1:nseqs
            uncompressed_oligos[m][i]=bases[m];
        end
    end
    return uncompressed_oligos
end

"""
    decompress_pool(compressed_pool)

Expand a pool of degenerate oligos into a single pool of non-degnerate oligos. 
"""
function decompress_pool(compressed_pool)
    decompressed_pool=decompress_oligo(compressed_pool[1])
    for i = 2:length(compressed_pool)
        subset=decompress_oligo(compressed_pool[i])
        for i in subset 
            push!(decompressed_pool,i)
        end 
    end 
    return decompressed_pool
end 

#Generate a random pool without replacement. Specify the randomer length, the size of the pool, and the available nucleotides.
#It can create ambiguous sequences, but it is possible that unique sequences will appear multiple times.

"""
    generate_random_pool(;randomerlen=6,randpoolsize=100,nucleotides=dna"ACGT")

Generate a random pool of oligos with length=randomerlen and pool size=randpoolsize using nucleotides
"""
function generate_random_pool(;randomerlen=6,randpoolsize=100,nucleotides=dna"ACGT")
    randpool=[]
    randomer=dna"-"^randomerlen
    for i = 1:randomerlen
        bases=shuffle(nucleotides)
        randomer[i]=bases[1];
    end
    push!(randpool,randomer)
    while length(randpool)<randpoolsize
        randomer=dna"-"^randomerlen
        for i = 1:randomerlen
            bases=shuffle(nucleotides)
            randomer[i]=bases[1];
        end
        check=true 
        for j = 1:length(randpool)
            if randomer==randpool[j]
                check=false
                break 
            end 
        end
        if check 
            push!(randpool,randomer)
        end
        if length(randpool)==length(nucleotides)^randomerlen
            break 
        end 
    end
    
    return randpool      
end 

# Given a pool of sequences all of length n, find all of the unambiguous sequences of length n that don't appear in that pool.

"""
    valid_candidates(sites)

Find all non-degenerate oligos that don't appear in sites. All oligos in sites must be of equal length.
"""
function valid_candidates(sites)
    lengths=length.(sites)
    if all(lengths[1].==lengths)
        n=length(sites[1])
        nsites=length(sites)
        candidates=generate_random_pool(;randomerlen=n,randpoolsize=4^n);
        for i=1:nsites
            candidates=candidates[findall(j -> !iscompatible(j,sites[i]),candidates)]
        end 
        return candidates
    else 
        println("provide an array of oligos of the same length")
        return nothing
    end
end  

# Given a list of valid candidates, determine the action space of those candidates. For each base position, return a list of compatible nucleotides. 
# If for example, the first postion does not contain an A or a T in any of the valid candidates, then the action space for that base position is dna"GCS". 
#returns bases, an array of dna sequences containing all of the valid actions for the agent to take at each position. 

"""
    define_action_space(candidates)

Find the available base codes for an oligo, such that only codes that appear in candidates at each postion are included.
"""
function define_action_space(candidates)
    if length(candidates)==0
        bases=[dna"" for i in 1:6]
        return bases
    else 
        randomerlen=length(candidates[1]);
        ncands=length(candidates)
        bases=[];
        for i=1:randomerlen 
            states=dna"-"^ncands;
            nucleotides=dna"----" # A list of nucleotides not present at a given position in the oligos, called "states" here 
            #println("TEST: $nucleotides")
            for j=1:ncands
                states[j]=candidates[j][i];
            end 
            if !occursin(dna"A",states)
                nucleotides[1]=DNA_A
            else 
                nucleotides[1]=DNA_Gap
            end 
            if !occursin(dna"T",states)
                nucleotides[2]=DNA_T
            else 
                nucleotides[2]=DNA_Gap
            end
            if !occursin(dna"C",states)
                nucleotides[3]=DNA_C
            else 
            nucleotides[3]=DNA_Gap
            end
            if !occursin(dna"G",states)
                nucleotides[4]=DNA_G
            else 
                nucleotides[4]=DNA_Gap
            end 
            #println("TEST: $states")
            push!(bases,omit_codes(nucleotides))
            #println("Test: $nucleotides")
        end
    end 
    return bases
end 

#given a list of base codes, return all of the base codes that do not contain any of those codes

""" 
    omit_codes(bases)

Helper function for define_action_space. 
"""
function omit_codes(bases::LongSequence{DNAAlphabet{4}})
    testcodes=dna"AGCTMRWSYKVHDBN"
    testcodes_copy=deepcopy(testcodes)
    indexarray=[]
    for i= 1:length(testcodes)
        count=0;
        for nucleotide in bases
            if iscompatible(nucleotide,testcodes[i])
                count+=1;
            end 
        end 
        if count >0
            push!(indexarray,i)
        end 
    end
    check=sort(indexarray,rev=true)
    for index in check
        deleteat!(testcodes_copy,index)   
    end 
   
    return testcodes_copy
end 

####################################################################
# Oligo Simulation Functions 
####################################################################


"""
    isvalid(seq, sites)

Check if any of the sites appear in the sequence. Also check for palindromes.

"""
function isvalid(seq, sites)
    if is_palindrome(seq)
        return false 
    end 
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
        candidates=shuffle(bases[i])
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

############################################################
#The follwoing functions are used as part of the various reward systems. They are called by each reward function when applicable  
###########################################################

function measure_hitcounts(randomer::LongSequence{DNAAlphabet{4}},mRNAs)
    mRNAnum=length(mRNAs);
    hitcounts=zeros(mRNAnum)
    for i=1:mRNAnum
        counter=0;
        indexcounter=0
        index=approxsearchindex(mRNAs[i][1:end],randomer,0); 
        indexcounter+=index;
        if index > 0
            counter +=1
        end 
        while index >0
            index=approxsearchindex(mRNAs[i][indexcounter+1:end],randomer,0);
            indexcounter+=index
            if index==0
                break
            else 
                counter+=1;
            end 
        end 
        hitcounts[i]=counter
    end 
    return hitcounts
end 


function mRNA_uniform_coverage_function_factory(randomer::LongSequence{DNAAlphabet{4}},mRNAs,hitcounts;uniformity_weight=1,total_hits_weight=1,kwargs...)
    hitcounts+=measure_hitcounts(randomer,mRNAs)
    function mRNA_uniform_coverage(randomer,mRNAs)
        mRNAnum=length(mRNAs);
        DULQ_score=0;
        new_hitcounts=measure_hitcounts(randomer,mRNAs)
        scored_hits=hitcounts+new_hitcounts;
        if sum(scored_hits)==0
            final_score=0
        else 
            sorted_hits=sort(scored_hits);
            Davg=mean(sorted_hits)
            lq=round(Integer,mRNAnum/4);
            Dlq=mean(sorted_hits[1:lq])
            DULQ_score=Dlq/Davg;
            filtered_hits=filter(x -> x > 0, scored_hits)
            final_score=length(filtered_hits)+total_hits_weight*sum(log10.(filtered_hits))+uniformity_weight*mRNAnum*DULQ_score
        end 
        return final_score
    end
    function mRNA_uniform_coverage_analysis(randomer,mRNAs)
        mRNAnum=length(mRNAs);
        DULQ_score=0;
        new_hitcounts=measure_hitcounts(randomer,mRNAs)
        scored_hits=hitcounts+new_hitcounts;
        if sum(scored_hits)==0
            final_score=0
        else 
            sorted_hits=sort(scored_hits);
            Davg=mean(sorted_hits)
            lq=round(Integer,mRNAnum/4);
            Dlq=mean(sorted_hits[1:lq])
            DULQ_score=Dlq/Davg;
            filtered_hits=filter(x -> x > 0, scored_hits)
            final_score=length(filtered_hits)+total_hits_weight*sum(log10.(filtered_hits))+uniformity_weight*mRNAnum*DULQ_score
        end 
        return final_score,hitcounts,DULQ_score, total_hits_weight, uniformity_weight
    end 
    return mRNA_uniform_coverage,mRNA_uniform_coverage_analysis,hitcounts
end 

function measure_hit_positions(randomer::LongSequence{DNAAlphabet{4}},mRNAs)
    mRNAnum=length(mRNAs)
    hit_positions=[[] for _ =1:mRNAnum]
    for i=1:mRNAnum
        indexcounter=0
        index=approxsearchindex(mRNAs[i][1:end],randomer,0); 
        indexcounter+=index;
        if index > 0
            push!(hit_positions[i],indexcounter)
        end 
        while index >0
            index=approxsearchindex(mRNAs[i][indexcounter+1:end],randomer,0);
            indexcounter+=index
            if index==0
                break
            else 
                push!(hit_positions[i],indexcounter)
            end 
        end 
    end 
    return hit_positions
end 

function calculate_gap_distances(hit_indicies,genelen,randomerlen)
    dists=[]
    if length(hit_indicies)>0  
        hit_indicies=sort(hit_indicies)
        push!(dists,hit_indicies[1])
        if length(hit_indicies)>1
            for i=2:length(hit_indicies)
                push!(dists,(hit_indicies[i]-hit_indicies[i-1]))
            end
        end
        push!(dists,(genelen-hit_indicies[end]))
    end 
    return dists
end 

#hit_positions= [[] for _ in 1:mRNAnum] 
function mRNA_within_gene_uniformity_function_factory(randomer::LongSequence{DNAAlphabet{4}},mRNAs,hit_positions;uniformity_weight=1,total_hits_weight=1,kwargs...)
    hit_positions=vcat.(hit_positions,measure_hit_positions(randomer,mRNAs));
    function mRNA_within_gene_uniformity(randomer,mRNAs)
        mRNAnum=length(mRNAs)
        gene_lengths=length.(mRNAs)
        new_hit_positions=measure_hit_positions(randomer,mRNAs)
        scored_hit_positions=vcat.(hit_positions,new_hit_positions)
        randomerlen=length(randomer)
        if sum(length.(scored_hit_positions))==0
            final_score=0
        else 
            gap_distances=calculate_gap_distances.(scored_hit_positions,gene_lengths,randomerlen)
            uniformity_score=zeros(mRNAnum)# "gap_standard_deviations_normalized" normalized by gene length
            for i =1:mRNAnum
                if length(gap_distances[i])==0
                    uniformity_score[i]=0
                else
                    sorted_gap_distances=sort(gap_distances);
                    Davg=mean(sorted_gap_distances)
                    lq=round(Integer,length(gap_distances)/4);
                    Dlq=mean(sorted_gap_distances[1:lq])
                    uniformity_score[i]=Dlq/Davg;
                end 
            end 
            filtered_hits=filter(x -> length(x)>0,scored_hit_positions)
            num_hits=length.(filtered_hits)
            final_score=length(filtered_hits)+ total_hits_weight*sum(log10.(num_hits))+ uniformity_weight*sum(uniformity_score)
        end 
        return final_score
    end 
    function mRNA_within_gene_uniformity_analysis(randomer,mRNAs)
        mRNAnum=length(mRNAs)
        gene_lengths=length.(mRNAs)
        new_hit_positions=measure_hit_positions(randomer,mRNAs)
        scored_hit_positions=vcat.(hit_positions,new_hit_positions)
        randomerlen=length(randomer)
        uniformity_score=zeros(mRNAnum)
        if sum(length.(scored_hit_positions))==0
            final_score=0
        else 
            gap_distances=calculate_gap_distances.(scored_hit_positions,gene_lengths,randomerlen)
            # "gap_standard_deviations_normalized" normalized by gene length
            for i =1:mRNAnum
                if length(gap_distances[i])==0
                    uniformity_score[i]=0
                else 
                    sorted_gap_distances=sort(gap_distances);
                    Davg=mean(sorted_gap+distances)
                    lq=round(Integer,length(gap_distances)/4);
                    Dlq=mean(sorted_gap_distances[1:lq])
                    uniformity_score[i]=Dlq/Davg;
                end 
            end 
            filtered_hits=filter(x -> length(x)>0,scored_hit_positions)
            num_hits=length.(filtered_hits)
            final_score=length(filtered_hits)+ total_hits_weight*sum(log10.(num_hits))+ uniformity_weight*sum(uniformity_score)
        end 
        return final_score,hit_positions,mean(uniformity_score),total_hits_weight,uniformity_weight
    end 
    return mRNA_within_gene_uniformity,mRNA_within_gene_uniformity_analysis, hit_positions
end 

function calculate_DULQ_score(array)
    if length(array)==0
        DULQ_score=0
    else 
        sorted_array=sort(array)
        grand_mean=mean(sorted_array)
        lq_idx=max(round(Integer,length(sorted_array)/4),1);
        lq_mean=mean(sorted_array[1:lq_idx])
        DULQ_score=lq_mean/grand_mean;
    end 
    return DULQ_score
end 
"""
    mRNA_all_rewards_function_factory(randomer::LongSequence{DNAAlphabet{4}},mRNAs,hitcounts,hit_positions;genes_hit_weight=1,intrauniformity_weight=0,total_hits_weight=0,interuniformity_weight=0,kwargs...)

Create reward and analysis functions to be used during rollout simulation.    

# Arguments
- `randomer`: Oligo to be added to the pool. The previously chosen cadidate.
- `mRNAs`: Array of mRNA sequences to be hit by NSR primers
- `hitcounts`: For each mRNA, the number of hits from previous oligos. 
- `hit_positions`: For each mRNA, the index of each hit from previous oligos
- `genes_hit_weight=1`: Weight for hitting every gene at least once 
- `total_hits_weight=0`: Weight for maximizing total hits
- `interuniformity_weight=0`: Weight for distributing equal numbers of hits to each gene
- `intrauniformity_weight=0`: Weight for placing hits uniformly along each gene. If 0, we do not track hit positions to reduce unnecessary computation.
"""
function mRNA_all_rewards_function_factory(randomer::LongSequence{DNAAlphabet{4}},mRNAs,hitcounts,hit_positions;genes_hit_weight=1,intrauniformity_weight=0,total_hits_weight=0,interuniformity_weight=0,kwargs...)
    if intrauniformity_weight==0  # we will only track the hit positions if the intrauniformity score matters. This will save unnecessary computation
        hitcounts+=measure_hitcounts(randomer,mRNAs)

        function mRNA_all_rewards_A(randomer,mRNAs)
            mRNAnum=length(mRNAs);
            interuniformity_score=0;
            new_hitcounts=measure_hitcounts(randomer,mRNAs)
            scored_hits=hitcounts+new_hitcounts;
            if sum(scored_hits)==0
                final_score=0
            else 
                filtered_hits=filter(x -> x > 0, scored_hits)
                interuniformity_score=calculate_DULQ_score(filtered_hits)
                final_score=genes_hit_weight*length(filtered_hits)+total_hits_weight*sum(filtered_hits)+interuniformity_weight*mRNAnum*interuniformity_score
            end 
            return final_score
        end
        function mRNA_all_rewards_analysis_A(randomer,mRNAs)
            mRNAnum=length(mRNAs);
            interuniformity_score=0;
            new_hitcounts=measure_hitcounts(randomer,mRNAs)
            scored_hits=hitcounts+new_hitcounts;
            genes_hit=0
            total_hits=0
            interuniformity_score=0
            intrauniformity_score=0
            if sum(scored_hits)==0
                final_score=0
                genes_hit=0
                total_hits=0
                interuniformity_score=0
                intrauniformity_score=0;
            else 
                filtered_hits=filter(x -> x > 0, scored_hits)
                interuniformity_score=calculate_DULQ_score(filtered_hits)
                genes_hit=length(filtered_hits)
                total_hits=sum(filtered_hits)
                intrauniformity_score= -1
                final_score=genes_hit_weight*length(filtered_hits)+total_hits_weight*sum(filtered_hits)+interuniformity_weight*mRNAnum*interuniformity_score
            end 
            return final_score,hitcounts,hit_positions,genes_hit,total_hits,interuniformity_score,mean(intrauniformity_score),genes_hit_weight,total_hits_weight,interuniformity_weight,intrauniformity_weight
        end 
        mRNA_all_rewards=mRNA_all_rewards_A
        mRNA_all_rewards_analysis=mRNA_all_rewards_analysis_A

    else ## intrauniformity_weight >0, keep track of hit positions
        hit_positions=vcat.(hit_positions,measure_hit_positions(randomer,mRNAs));
        hitcounts=length.(hit_positions)
        function mRNA_all_rewards_B(randomer,mRNAs)
            mRNAnum=length(mRNAs)
            gene_lengths=length.(mRNAs)
            new_hit_positions=measure_hit_positions(randomer,mRNAs)
            scored_hit_positions=vcat.(hit_positions,new_hit_positions)
            randomerlen=length(randomer)
            if sum(length.(scored_hit_positions))==0
                final_score=0
            else 
                gap_distances=calculate_gap_distances.(scored_hit_positions,gene_lengths,randomerlen)
                gap_distances=filter(x -> length(x)>0,gap_distances)
                intrauniformity_scores=calculate_DULQ_score.(gap_distances)
                filtered_hit_positions=filter(x -> length(x)>0,scored_hit_positions)
                filtered_hitcounts=(length.(filtered_hit_positions));
                interuniformity_score=calculate_DULQ_score(filtered_hitcounts)
                filtered_hits=length.(filtered_hit_positions)
                final_score=genes_hit_weight*length(filtered_hits)+ total_hits_weight*sum(filtered_hits)+ interuniformity_weight*mRNAnum*interuniformity_score+intrauniformity_weight*mRNAnum*mean(intrauniformity_scores)
            end 
            return final_score
        end 
        function mRNA_all_rewards_analysis_B(randomer,mRNAs)
            mRNAnum=length(mRNAs)
            gene_lengths=length.(mRNAs)
            new_hit_positions=measure_hit_positions(randomer,mRNAs)
            scored_hit_positions=vcat.(hit_positions,new_hit_positions)
            randomerlen=length(randomer)
            genes_hit=0
            total_hits=0
            interuniformity_score=0
            intrauniformity_scores=0
            if sum(length.(scored_hit_positions))==0
                final_score=0
                genes_hit=0
                total_hits=0
                interuniformity_score=0
                intrauniformity_scores=0
            else 
                gap_distances=calculate_gap_distances.(scored_hit_positions,gene_lengths,randomerlen)
                gap_distances=filter(x -> length(x)>0,gap_distances)
                intrauniformity_scores=calculate_DULQ_score.(gap_distances)
                filtered_hit_positions=filter(x -> length(x)>0,scored_hit_positions)
                filtered_hitcounts=(length.(filtered_hit_positions));
                interuniformity_score=calculate_DULQ_score(filtered_hitcounts)
                filtered_hits=length.(filtered_hit_positions)
                genes_hit=length(filtered_hits)
                total_hits=sum(filtered_hits)
                final_score=genes_hit_weight*length(filtered_hits)+ total_hits_weight*sum(filtered_hits)+ interuniformity_weight*mRNAnum*interuniformity_score+intrauniformity_weight*mRNAnum*mean(intrauniformity_scores)
            end 
            return final_score,hitcounts,hit_positions,genes_hit,total_hits,interuniformity_score,mean(intrauniformity_scores),genes_hit_weight,total_hits_weight,interuniformity_weight,intrauniformity_weight
        end
        mRNA_all_rewards=mRNA_all_rewards_B
        mRNA_all_rewards_analysis=mRNA_all_rewards_analysis_B
    end 
    return mRNA_all_rewards,mRNA_all_rewards_analysis,hitcounts,hit_positions 
end 
#################################################
# Simulation function. Call Reward functions and simulate a valid randomer 
##################################################
function simulate(prefix,bases,sites,mRNAs,reward_function; nsims=1000, horizon=length(bases), kwargs...)
    rewards=zeros(nsims)
    true_horizon = min(horizon, length(bases))
    n=length(prefix)+true_horizon;
    sim_randomer_list=[dna"-"^n for _ in 1:nsims];
    
    Threads.@threads for i = 1:nsims   #This is where multithreading is incorporated 
        sim_randomer_list[i]=find_valid_randomer(prefix,bases[1:true_horizon],sites);
        rewards[i]= reward_function(sim_randomer_list[i],mRNAs)
    end
    max_sim_score,max_index=findmax(rewards);
    max_sim_randomer=sim_randomer_list[max_index];
    return mean(rewards),max_sim_randomer,max_sim_score

end

########################################
# Rollout Function
#######################################
function oligo_rollout(bases,sites, mRNA; reward_function=mRNA_uniform_coverage, kwargs...)
    n=length(bases)
    randomer= dna"-" ^ n
    # Find the longest blocking site and use this as the horizon. (they should all be equal)
    max_len = map(length, sites) |> maximum
    candidates=valid_candidates(sites);
    best_sim_randomer=dna"-"^n
    best_sim_score=0;
    for i= 1:n
        best_base= DNA_Gap
        best_score=0;
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        bases=define_action_space(candidates[findall(j -> iscompatible(j[start:i-1],randomer[start:i-1]),candidates)])
        if length(bases[i])==1
            best_base=bases[i][1];
        else 
            for base in bases[i]
                if ~isvalid(randomer[start:i-1] * base, sites)
                    # adding this bases creates a restriction site; skip it
                    continue
                end
                score,max_sim_randomer,max_sim_score= simulate(randomer[start:i-1] * base, bases[i+1:end], sites, mRNA,reward_function;kwargs...)
                # Update use this base if:
                #   1. the mean randomer mRNA coverage is higher than the previous best, OR
                #   2. the randomer mRNA coverage is tied with the previous best, but the new base has
                #      higher individual degeneracy           
                if score > best_score || (score == best_score && degeneracy(base) > degeneracy(best_base))
                    best_score = score
                    best_base = base
                end
                ## Update the overall best simulation result if 
                #   1. the best globally simulated randomer mRNA score is better than the previous best
                #   2. the best globally simulated randomer mRNA score is tied but the randomer has higher degeneracy
                if max_sim_score >best_sim_score || (max_sim_score==best_sim_score && degeneracy(max_sim_randomer)>degeneracy(best_sim_randomer))
                    best_sim_score=max_sim_score;
                    best_sim_randomer=max_sim_randomer
                end
            end
        end 
        if best_base == DNA_Gap
            # The sequence has terminated; there are no bases that can be added without making
            # a restriction site.
            break
        else
            randomer[i] = best_base
            
        end
        #println("randomer: $randomer")
    end
    final_randomer_score=reward_function(randomer,mRNA)
    if best_sim_score >final_randomer_score
        randomer=best_sim_randomer
        println("Triggered Best Simulated Randomer")
    end 
    return randomer
end

####################################
#Functions for Reading Genome Data CSV Files Imported from R
######################################

## Read in the list of sequences you want to exclude in the designed primers. 
# Function finds every unique hexamer in the list of sequences and stores it as an array
#. For example, use list of rRNA and tRNA sequences. This funciton takes care of orientation. We end up using the reverse complement since we are finding primers that won't hybridize 
## Example
#blocking_sites=read_blocking_sites(CSV.read("./SMU_UA159_rRNA_tRNA.csv"))
function read_blocking_sites(data; randomerlen=6)
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
# RNAs=read_rRNA_tRNA_list(CSV.read("./SMU_UA159_rRNA_tRNA.csv"))
#This function Does not hash up the rRNAs and tRNAs like read_blocking_sites does. That is the only difference
function read_rRNA_tRNA_list(data)
    RNAs=[]
    for i=1:nrow(data)
        seq=reverse_complement(LongDNASeq(string(data[i,:Sequence])));
        if isvalid(seq,RNAs)
            push!(RNAs,seq)
        end 
    end 
    return RNAs
end 
## Read in the list of sequences you want the designed primers to hybridize to
# Function takes care of the orientation, we want to end up using the reverse complmement. 
## Example
# mRNAs=read_mRNA_list(CSV.read("./SMU_UA159_Genes.csv"))
function read_mRNA_list(data)
    mRNAs=[]
    for i=1:nrow(data)
        seq=reverse_complement(LongDNASeq(string(data[i,:mRNA_seq])));
        if isvalid(seq,mRNAs)
            push!(mRNAs,seq)
        end 
    end 
    return mRNAs
end

#Finds the number of bases as well as the GC content of the given file 
# Example 
# bases,GC_content=calculate_transcriptome_stats(read_mRNA_list(CSV.read("./SMU_UA159_Genes.csv")))
function calculate_transcriptome_stats(sequences)
    bases=0;  #number of bases in the file 
    GC_count=0 # Count of G's and C's in file
    for i=1:length(sequences)
        bases += length(sequences[i]);
        comp=composition(sequences[i]);
        GC_count += comp[DNA_G] +comp[DNA_C];
    end 
    GC_content= GC_count/bases*100;
    return bases,GC_content
end 





#################################################
#### Extra Analysis Tools 
#################################################

## Finds the edit distance distribution for a given designed primer to all of the blocking sites 
function base_error_dist(randomer,sites)
    error_dist=[]
    for i in 0:length(randomer)
        count=0
        found_sites=[]
        for j in 1:length(sites)
            valid=approxsearchindex(randomer,sites[j],i)
            if valid == 1
                count += valid
                push!(found_sites,j)
            end
        end
        deleteat!(sites,found_sites)
        push!(error_dist,count)
    end

    return error_dist
end

function ID_genes_missing(data,genes_missing)
    missing_ID=[];
    for i in genes_missing
        push!(missing_ID,data[i,:New_ID])
    end 
    return missing_ID
end 


## Calculates the distribution number of randomers that hit each gene at least once. Not very useful for most cases 
function coverage_distribution(randomers, mRNA)
    coverage_dist=zeros(length(randomers)+1)
     for seq in mRNA
         count=0;
         for randomer in randomers
             if occursin(randomer, seq)
                 count+=1;
             end
         end
         coverage_dist[count+1]+=1;
 
     end
     return coverage_dist
end

# Calcualtes the number of times each gene is hit by all randomers in the pool. Useful for comparing the "between gene uniformity". 
function gene_coverage_counts(randomers,mRNAs)
    gene_hits=zeros(length(mRNAs))
    randomerlen=length(randomers[1]);
    for i in 1:length(mRNAs)
        counter=0;
        for j in 1:length(mRNAs[i])-randomerlen
            for k in 1:length(randomers)
                if occursin(randomers[k],mRNAs[i][j:j+randomerlen-1])
                    counter+=1
                end
            end
        end
        gene_hits[i]=counter
    end
    data=DataFrame(counts=gene_hits)
    return data

end



 ########################################################################
 # MAIN RUNNING FUNCTION
 ########################################################################

 #Quick Test#
 # data=run_NSR_RL(pool_size=3,nsims=50)
 #CSV.write("NSR_RL_S_mutans_4_12_21",data)

 #Input the Species name, the length of primer, and the size of the pool. 
function run_NSR_RL(;species="S_mutans",randomerlen=6,pool_size=25,all_bases = dna"AGCTMRWSYKVHDBN",kwargs...)
    
    #Intialization
    randomers=[]
    mRNAfile=string("./Genome Data/",species,"_Genes.csv")
    rRNAfile=string("./Genome Data/",species,"_rRNA_tRNA.csv")
    mRNA_DF=CSV.read(mRNAfile,DataFrame);
    rRNA_DF=CSV.read(rRNAfile,DataFrame);
    mRNAs=read_mRNA_list(mRNA_DF);
    rRNAs=read_rRNA_tRNA_list(rRNA_DF);
    numgenes=length(mRNAs)
    blocking_sites=read_blocking_sites(rRNA_DF,randomerlen=randomerlen);
    mRNAbp,mRNA_GC=calculate_transcriptome_stats(mRNAs);
    rRNAbp,rRNA_GC=calculate_transcriptome_stats(rRNAs);
    transcriptome_size=mRNAbp+rRNAbp;
    all_ns = [all_bases for _ in 1:randomerlen];
    times=[]
    num_blocked_sites=zeros(pool_size)
    println("Initialization Successful")
    ##Initialize the reward-associated data structure and reward function here. ex hitcounts=zeros(numgenes), hit_positions=[[] for _ in 1:numgenes]
    hitcounts=zeros(numgenes)
    hit_positions=[[] for _ in 1:numgenes]
    #Call the appropriate reward function factory 
    reward_function,analysis_function,hitcounts,hit_positions=mRNA_all_rewards_function_factory(dna"-"^randomerlen,mRNAs,hitcounts,hit_positions;kwargs...)
    ##
    for i = 1:pool_size
        num_blocked_sites[i]=length(blocking_sites)
        randomerstats= @timed oligo_rollout(all_ns,blocking_sites,mRNAs; reward_function=reward_function,kwargs...); #create new randomer and add to pool
        push!(randomers,randomerstats.value)
        push!(times,randomerstats.time)  
        push!(blocking_sites,LongDNASeq(randomers[i]));
        push!(blocking_sites,LongDNASeq(reverse_complement(randomers[i]))); ## avoid adding randomers that are reverse complements of each other
        blocking_sites=unique(decompress_pool(blocking_sites));
        reward_function,analysis_function,hitcounts,hit_positions=mRNA_all_rewards_function_factory(randomers[i],mRNAs,hitcounts,hit_positions; kwargs...) #officially add hit positions from the new randomer to the reward function 
        #println("Iteration $i")
    end 
    final_score,hitcounts,hit_positions,genes_hit,total_hits,interuniformity_score,intrauniformity_score,genes_hit_weight,total_hits_weight,interuniformity_weight,intrauniformity_weight =analysis_function(dna"-"^randomerlen,mRNAs)
    genes_missing=findall(x -> x==0,hitcounts)
    missing_IDs=join(ID_genes_missing(mRNA_DF,genes_missing),",")
    totaldeg=pooled_degeneracy(randomers)
    total_time=sum(times)
    string_time=join(times,",")
    string_num_blocked_sites=join(num_blocked_sites,",")
    string_randomers=join(convert.(String,randomers),",")
    data=DataFrame(poolsize=pool_size, primer_length=randomerlen, genes_hit_weight=genes_hit_weight, total_hits_weight=total_hits_weight,interuniformity_weight=interuniformity_weight,intrauniformity_weight=intrauniformity_weight,species=species,num_genes=numgenes,
        transcriptome_size=transcriptome_size,mRNA_GC=mRNA_GC,blocked_GC=rRNA_GC,randomers=string_randomers,degeneracy=totaldeg,num_genes_hit=genes_hit,
        genes_missing=missing_IDs,total_hits=total_hits,interuniformity_score=interuniformity_score,intrauniformity_score=intrauniformity_score,individual_times=string_time,time=total_time,num_blocked_sites=string_num_blocked_sites,final_score=final_score)
    return data 
end 

function benchmark_NSR_RL()
    DF1=DataFrame()
    nexps=1;
    genes_hit_weights=[1]
    total_hits_weights=[1e-4]
    interuniformity_weights=[1]
    intrauniformity_weights=[1]
    species="S_mutans"
    filename="./S_mutans_NSR_Pool_7_22_21.csv"
    for i = 1:nexps
        data=run_NSR_RL(;species=species,pool_size=100,nsims=100,genes_hit_weight=genes_hit_weights[i],total_hits_weight=total_hits_weights[i],interuniformity_weight=interuniformity_weights[i],intrauniformity_weight=intrauniformity_weights[i])
        DF1=vcat(DF1,data);
        CSV.write(filename,DF1)
    end 
end


##Code for analyzing previously made primer pools. stores the same data that run_NSR_RL stores 
function analyze_NSR_Pool(randomers;species="S_mutans",mRNAfile="./SMU_UA159_Genes.csv",rRNAfile="./SMU_UA159_rRNA_tRNA.csv",kwargs...)
    mRNAfile=string("./Genome Data/",species,"_Genes.csv")
    rRNAfile=string("./Genome Data/",species,"_rRNA_tRNA.csv")
    mRNA_DF=CSV.read(mRNAfile,DataFrame);
    rRNA_DF=CSV.read(rRNAfile,DataFrame);
    mRNAs=read_mRNA_list(mRNA_DF);
    rRNAs=read_rRNA_tRNA_list(rRNA_DF);
    numgenes=length(mRNAs)
    blocking_sites=read_blocking_sites(rRNA_DF,randomerlen=randomerlen);
    randomerlen=length(randomers[1]);
    mRNAbp,mRNA_GC=calculate_transcriptome_stats(mRNAs);
    rRNAbp,rRNA_GC=calculate_transcriptome_stats(rRNAs);
    transcriptome_size=mRNAbp+rRNAbp;
    hitcounts=zeros(numgenes)
    hit_positions=[[] for _ in 1:numgenes]
    reward_function,analysis_function,hitcounts,hit_positions=mRNA_all_rewards_function_factory(dna"-"^randomerlen,mRNAs,hitcounts,hit_positions;genes_hit_weight=1,total_hits_weight=1,interuniformity_weight=1,intrauniformity_weight=1)
    for i=1:length(randomers)
        reward_function,analysis_function,hitcounts,hit_positions=mRNA_all_rewards_function_factory(randomers[i],mRNAs,hitcounts,hit_positions;genes_hit_weight=1,total_hits_weight=1,interuniformity_weight=1,intrauniformity_weight=1); #officially add hit positions from the new randomer 
    end
    final_score,hitcounts,hit_positions,genes_hit,total_hits,interuniformity_score,intrauniformity_score,genes_hit_weight,total_hits_weight,interuniformity_weight,intrauniformity_weight =analysis_function(dna"-"^randomerlen,mRNAs)
    genes_missing=findall(x->x==0,hitcounts);
    missing_IDs=join(ID_genes_missing(mRNAdata,genes_missing),",")
    totaldeg=pooled_degeneracy(randomers)
    total_time=0
    string_time="N/A"
    string_randomers=join(convert.(String,randomers),",")
    genes_hit_weight=0
    total_hits_weight=0
    interuniformity_weight=0
    intrauniformity_weight=0
    pool_size=length(randomers);
    data=DataFrame(poolsize=pool_size, primer_length=randomerlen, genes_hit_weight=genes_hit_weight, total_hits_weight=total_hits_weight,interuniformity_weight=interuniformity_weight,intrauniformity_weight=intrauniformity_weight,species=species,num_genes=numgenes,
        transcriptome_size=transcriptome_size,mRNA_GC=mRNA_GC,blocked_GC=rRNA_GC,randomers=string_randomers,degeneracy=totaldeg,num_genes_hit=genes_hit,
        genes_missing=missing_IDs,total_hits=total_hits,interuniformity_score=interuniformity_score,intrauniformity_score=intrauniformity_score,individual_times=string_time,time=total_time)
    return data
end 

function analyze_all_genomes_NSR_Pools(species_list)
    DF1=DataFrame()
    n=length(species_list)
    for i = 1:n
        randomers=LongDNASeq.(CSV.read(string("./Brute Force Pools V2 Compressed/",species_list[i],"_NSR_Brute_Force_Compressed.csv"),DataFrame)[:Sequence])
        mRNAfile=string("./Genome Data/",species_list[i],"_Genes.csv")
        rRNAfile=string("./Genome Data/",species_list[i],"_rRNA_tRNA.csv")
        data=analyze_NSR_Pool(randomers;species=species_list[i],mRNAfile=mRNAfile,rRNAfile=rRNAfile)
        DF1=vcat(DF1,data);
        pct_complete=i/n*100;
        println("$pct_complete % Complete")
        CSV.write("./Brute_Force_Pool_Compressed_Analysis_all_rewards.csv",DF1)
    end 
    
    
end 
function analyze_NSR_Pool_gene_coverage(species_name)
    randomers=LongDNASeq.(CSV.read(string("./Brute Force Pools V2/",species_name,"_NSR_Brute_Force_all_randomers.csv"),DataFrame)[:Sequence])
    mRNAfile=string("./Genome Data/",species_name,"_Genes.csv")
    data=gene_coverage_counts(randomers,mRNAfile);
    outputfile=string("./Brute Force Pools V2 Gene Coverage/",species_name,"_Gene_Coverage.csv")
    CSV.write(outputfile,data)
end 

function analyze_cumulative_NSR_Pool(randomers,species_name;outputfile="./S_mutans_Cumulative_NSR_Pool.csv" ,kwargs...)
    mRNAfile=string("./Genome Data/",species_name,"_Genes.csv")
    mRNAdata=CSV.read(mRNAfile,DataFrame);
    mRNAs=read_mRNA_list(mRNAdata);
    nrandomers=length(randomers);
    randomerlen=length(randomers[1]);
    cumulative_hits=zeros(nrandomers)
    cumulative_score=zeros(nrandomers)
    cumulative_interuniformity_score=zeros(nrandomers)
    cumulative_intrauniformity_score=zeros(nrandomers)
    randomer_num=1:nrandomers;
    str_randomers=convert.(String,randomers);
    cumulative_genes_hit=zeros(nrandomers);
    hitcounts=zeros(length(mRNAs))
    hit_positions=[[] for _ in 1:length(mRNAs)]
    reward_function,analysis_function,hitcounts,hit_positions=mRNA_all_rewards_function_factory(dna"-"^randomerlen,mRNAs,hitcounts,hit_positions;genes_hit_weight=1,total_hits_weight=1,interuniformity_weight=1,intrauniformity_weight=1)
    for i=1:nrandomers
        reward_function,analysis_function,hitcounts,hit_positions=mRNA_all_rewards_function_factory(randomers[i],mRNAs,hitcounts,hit_positions;genes_hit_weight=1,total_hits_weight=1,interuniformity_weight=1,intrauniformity_weight=1); #officially add hit positions from the new randomer 
        final_score,hitcounts,hit_positions,genes_hit,total_hits,interuniformity_score,intrauniformity_score,genes_hit_weight,total_hits_weight,interuniformity_weight,intrauniformity_weight=analysis_function(dna"-"^randomerlen,mRNAs)
        cumulative_score[i]=final_score
        cumulative_hits[i]=total_hits
        cumulative_interuniformity_score[i]=interuniformity_score
        cumulative_genes_hit[i]=genes_hit
        cumulative_intrauniformity_score[i]=intrauniformity_score
    end 
    data=DataFrame(run=randomer_num,randomer=str_randomers,cumulative_genes_hit=cumulative_genes_hit,cumulative_hits=cumulative_hits,cumulative_interuniformity_score=cumulative_interuniformity_score,cumulative_intrauniformity_score=cumulative_intrauniformity_score,cumulative_score=cumulative_score)
    CSV.write(outputfile,data)
end 

function cumulative_runtime()
    data=CSV.read("./S_mutans_NSR_Refactor_NSR_RL_All_Rewards_BF_Comparison_6_21_21.csv",DataFrame)
    hitcountstime=parse.(Float64,(split(data[2,19],",")))
    hit_pos_time=parse.(Float64,(split(data[4,19],",")))
    run=1:length(hitcountstime)
    cum_hitcountstime=zeros(length(hitcountstime))
    cum_hit_pos_time=zeros(length(hit_pos_time))
    for i =1:length(hit_pos_time)
        cum_hitcountstime[i]=sum(hitcountstime[1:i])
        cum_hit_pos_time[i]=sum(hit_pos_time[1:i])
    end 
    df=DataFrame(run=run,without_intra=cum_hitcountstime,with_intra=cum_hit_pos_time)
    CSV.write("NSR_pool_size_timing_6_23_21.csv",df)
    end 



#data=CSV.read("NSR_RL/Brute Force Pools V2 Compressed/S_mutans_NSR_Brute_Force_Compressed.csv",DataFrame)
#randomers=LongDNASeq.(data[:Sequence])
#analyze_cumulative_NSR_Pool(randomers,"S_mutans";outputfile="./NSR_RL/Experiments/S_mutans_Cumulative_BF_Compressed.csv")


NSR_RL/Brute Force Pools V2 Compressed/S_mutans_NSR_Brute_Force_Compressed.csv