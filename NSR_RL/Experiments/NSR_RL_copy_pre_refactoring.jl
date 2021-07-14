import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences, FASTX, XLSX, CPUTime

#########
#The following functions are the base set of tools that this program uses 

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

Check if any of the sites appear in the sequence. Also check for palindromes by calling is palindrom

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
function iscompatible_sequence(seqa,seqb)
    n=length(seqa)
    count=0
    for i=1:n
        if iscompatible(seqa[i],seqb[i])
            count+=1;
        end
    end 
    return count==n
end 

function decompress_oligo(comp_oligo)
    basedict=Dict(DNA_A => dna"A", DNA_T => dna"T", DNA_C => dna"C", DNA_G => dna"G",DNA_R =>dna"AG",DNA_Y=>dna"CT",DNA_S=>dna"GC",
    DNA_W=>dna"AT",DNA_K=>dna"GT",DNA_M=>dna"AC",DNA_B=>dna"CGT",DNA_D=>dna"AGT",DNA_H=>dna"ACT",DNA_V=>dna"ACG",DNA_N=>dna"ACGT",DNA_Gap=>dna"-")
    nseqs=trunc(Int128,degeneracy(comp_oligo;uselog2=false))
    randomerlen=length(comp_oligo)
    blockarray=zeros(Int128,randomerlen+1)
    blockarray[1]=nseqs
    degsarray=zeros(Int128,randomerlen)
    repsarray=zeros(Int128,randomerlen+1)
    repsarray[1]=1
    uncomp_oligos=[dna"-"^randomerlen for _ in 1:nseqs]
    for i = 1:randomerlen 
        nucleotides=basedict[comp_oligo[i]]
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
            uncomp_oligos[m][i]=bases[m];
        end
    end
    return uncomp_oligos
end
function decompress_pool(compressed_pool)
    uncomppool=decompress_oligo(compressed_pool[1])
    for i = 2:length(compressed_pool)
        subset=decompress_oligo(compressed_pool[i])
        for i in subset 
            push!(uncomppool,i)
        end 
    end 
    return uncomppool 
end 

#####################
#Define the action space 
#####################
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

function valid_candidates(sites)
    n=length(sites[1])
    nsites=length(sites)
    candidates=generate_random_pool(;randomerlen=n,randpoolsize=4^n);
    for i=1:nsites
        candidates=candidates[findall(j -> !iscompatible_sequence(j,sites[i]),candidates)]
    end 
    return candidates
end 


function define_action_space(sites)
    if length(sites)==0
        bases=[dna"" for i in 1:6]
        return bases
    else 
        randomerlen=length(sites[1]);
        nsites=length(sites)
        bases=[];
        for i=1:randomerlen 
            states=dna"-"^nsites;
            nucleotides=dna"----"
            #println("TEST: $nucleotides")
            for j=1:nsites
                states[j]=sites[j][i];
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
function omit_codes(bases)
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




#######################
#The follwoing functions are used as part of the various reward systems 
######################


#SCORING SYSTEM FOR HITTING ALL mRNAs. Not used in the weighted scoring system 
function mRNA_coverage(randomer,mRNA)
    hits=0;
    for seq in mRNA 
        if occursin(randomer,seq)
            hits+=1;
        end 
    end
    return hits
end

#Weighted Scoring System for mRNA coverage
function mRNA_weighted(randomer,mRNA,weights)
score=0;
for i=1:length(mRNA)
    if occursin(randomer,mRNA[i])
        score+=weights[i];
    end
end
return score
end

# Keeps track of the number of times a given gene has been hit. It is used to compute the weights above 
function mRNA_hitcounter(randomer,mRNA,hitcounts)
    for i=1:length(mRNA)
        if occursin(randomer,mRNA[i])
            hitcounts[i]+=1;
        end
    end
    return hitcounts
end
function mRNA_uniform_coverage(randomer,mRNAs;hitcounts=zeros(length(mRNAs)),uniformity_weight=1,total_hits_weight=1,overwrite_hitcounts=false, kwargs...)
    mRNAnum=length(mRNAs);
    gene_scores=zeros(mRNAnum);
    dist_score=0;
    randomerlen=length(randomer);
    new_hitcounts=zeros(mRNAnum)
    #println("hitcounts: $hitcounts")
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
     

        
        new_hitcounts[i]=counter
        hit_score=hitcounts[i]+new_hitcounts[i]
        if hit_score==0
            gene_scores[i]=0;
        else 
            gene_scores[i]=1+total_hits_weight*log10(hit_score)
        end
    end 
    scored_hits=hitcounts+new_hitcounts;
    if sum(scored_hits)==0
        final_score=0
    else 
        sorted_hits=sort(scored_hits);
        DUavg=mean(sorted_hits)
        lq=round(Integer,mRNAnum/4);
        DUlq=mean(sorted_hits[1:lq])
        dist_score=DUlq/DUavg;
        final_score=sum(gene_scores)+uniformity_weight*mRNAnum*dist_score;
    end 
    if overwrite_hitcounts
        return final_score,scored_hits,dist_score, uniformity_weight,total_hits_weight
    else 
        return final_score,hitcounts,dist_score,uniformity_weight,total_hits_weight
    end 
end 




function mRNA_multiobjective(randomer, mRNAs; indicies=[[] for i=1:10], b1=1,b2=1,overwrite_indicies=false, return_indicies=false,kwargs...)
    # If overwrite_indicies=false, the funciton will score the randomer as if the randomer were to be added, but will not actually add the randomer
    # If overwrite_indicies=true, then the function will add the randomer hit sites to the offical count 
    # b1 = weight for uniformity score 
    # b2 = weight for hits per unit length 
    mRNAnum=length(mRNAs)
    gene_scores=zeros(mRNAnum);
    means=zeros(mRNAnum);
    stds=ones(mRNAnum);
    randomerlen=length(randomer);
    indicies_copy=deepcopy(indicies);
    for i=1:mRNAnum
        iteration=0
        index=approxsearchindex(mRNAs[i][1:end],randomer,0);  
        if index > 0
            push!(indicies_copy[i],index)
        end
        iteration = iteration+1;
        while index > 0
            index=approxsearchindex(mRNAs[i][((indicies_copy[i][end]+1):end)],randomer,0);
            iteration=iteration+1;
            if index == 0
                break
            else 
                push!(indicies_copy[i],index+indicies_copy[i][end])
                
            end
        end
        if length(indicies_copy[i])==0
            continue
        end 
        dists=calculate_gap_distances(indicies_copy[i],length(mRNAs[i]),randomerlen)
        stds[i]=std(dists)/length(mRNAs[i])
        gene_scores[i]=1+b1*(1-stds[i])+b2*length(indicies_copy[i])/length(mRNAs[i]);  
    end

    if overwrite_indicies && return_indicies
    return indicies_copy, b1, b2
    elseif ~overwrite_indicies &&  return_indicies
    return  indicies
    else  # by having return_indicies=false, you can just return the gene scores. This just makes the output of the function more flexible. 
        return sum(gene_scores)
    end
end

function calculate_gap_distances(hit_indicies,genelen,randomerlen)
    dists=[]
    hit_indicies=sort(hit_indicies)
    push!(dists,hit_indicies[1])
    if length(hit_indicies)>1
        for i=2:length(hit_indicies)
            push!(dists,(hit_indicies[i]-hit_indicies[i-1]-randomerlen))
        end
    end
    push!(dists,(genelen-hit_indicies[end]-randomerlen))
    return dists
end 

# Rollout reward systems
#3 factor scoring policy. This is the default reward system 
function simulate_random_mRNA_uniform_coverage(prefix,bases,sites,mRNAs; nsims=1000, horizon=length(bases), kwargs...)
    mRNA_score=zeros(nsims)
    true_horizon = min(horizon, length(bases))
    n=length(prefix)+true_horizon;
    sim_randomer_list=[dna"-"^n for _ in 1:nsims];
    
    Threads.@threads for i = 1:nsims   #This is where multithreading is incorporated 
        sim_randomer_list[i]=find_valid_randomer(prefix,bases[1:true_horizon],sites);
        mRNA_score[i],b,c,d,e= mRNA_uniform_coverage(sim_randomer_list[i],mRNAs;overwrite_hitcounts=false, kwargs...)
    end
    max_sim_score,max_index=findmax(mRNA_score);
    max_sim_randomer=sim_randomer_list[max_index];
    return mean(mRNA_score),max_sim_randomer,max_sim_score

end
function simulate_random_mRNA_multiobjective(prefix,bases,sites,mRNAs; nsims=1000, horizon=length(bases), kwargs...)
    mRNA_score=zeros(nsims)
    true_horizon = min(horizon, length(bases))
    n=length(prefix)+true_horizon;
    sim_randomer_list=[dna"-"^n for _ in 1:nsims];
    Threads.@threads for i = 1:nsims   #This is where multithreading is incorporated 
        sim_randomer_list[i]=find_valid_randomer(prefix,bases[1:true_horizon],sites);
        mRNA_score[i]= mRNA_multiobjective(sim_randomer_list[i],mRNAs;return_indicies=false, kwargs...)
    end
    max_sim_score,max_index=findmax(mRNA_score);
    max_sim_randomer=sim_randomer_list[max_index];
    return mean(mRNA_score),max_sim_randomer,max_sim_score
end

#Weighted scoring policy. The score for hitting a given gene decreases according to the number of times it's been hit
function simulate_random_mRNA_weighted(prefix,bases,sites,mRNA; nsims=1000, horizon=length(bases))
    mRNA_score=zeros(nsims)
    true_horizon = min(horizon, length(bases))
    Threads.@threads for i = 1:nsims
        mRNA_score[i]= mRNA_weighted(find_valid_randomer(prefix,bases[1:true_horizon],sites),mRNA,weights)
    end
    return mean(mRNA_score)
end

# Unweighted scoring policy. Finds randomers that hit the most genes at least once 
function simulate_random_mRNA(prefix,bases,sites,mRNA; nsims=1000, horizon=length(bases))
    mRNA_hits=zeros(nsims)
    true_horizon = min(horizon, length(bases))
    Threads.@threads for i = 1:nsims
        mRNA_hits[i]= mRNA_coverage(find_valid_randomer(prefix,bases[1:true_horizon],sites),mRNA)
    end
    return mean(mRNA_hits)
end


### MAIN FUNCTION
function rollout_rt(bases,sites, mRNA; simulate=simulate_random_mRNA_multiobjective, kwargs...)
    n=length(bases)
    randomer= dna"-" ^ n
    # Find the longest blocking site and use this as the horizon. (they should all be equal)
    max_len = map(length, sites) |> maximum
    horizon = max_len
    candidates=valid_candidates(sites);
    checkpoint=length(candidates)
    best_sim_randomer=dna"-"^n
    best_sim_score=0;
    #println("Valid Candidates: $checkpoint")
    for i= 1:n
        best_base= DNA_Gap
        best_score=0;
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        bases=define_action_space(candidates[findall(j -> iscompatible_sequence(j[start:i-1],randomer[start:i-1]),candidates)])
        #println("Bases $bases")
        if length(bases[i])==1
            best_base=bases[i][1];
        else 
            for base in bases[i]
                if ~isvalid(randomer[start:i-1] * base, sites)
                    # adding this bases creates a restriction site; skip it
                    continue
                end
                score,max_sim_randomer,max_sim_score= simulate(randomer[start:i-1] * base, bases[i+1:end], sites, mRNA;kwargs...)
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
    final_randomer_score,hitcounts,dist_score,uniformity_weight,total_hits_weight=mRNA_uniform_coverage(randomer,mRNA;kwargs...)
    if best_sim_score >final_randomer_score
        randomer=best_sim_randomer
        println("Triggered Best Simulated Randomer")
    end 
    return randomer
end


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






#### Extra Analysis Funcitons 

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
# genes_hit,genes_missing,total_hits,uniformity_score=hit_stats(indicies,mRNAs,6)
function mRNA_multiobjective_hit_stats(indicies,mRNAs,randomerlen)
    genes_hit=0;
    genes_missing=[]
    total_hits=0;
    stds=ones(length(indicies))
    for i=1:length(indicies);
        if length(indicies[i])>0
            genes_hit+=1;
            dists=calculate_gap_distances(indicies[i],length(mRNAs[i]),randomerlen)
            stds[i]=std(dists)/length(mRNAs[i]);
        else 
            push!(genes_missing,i)
        end 
        total_hits+=length(indicies[i]);
    end 
    return genes_hit,genes_missing,total_hits,mean(1 .- stds)
end 


function ID_genes_missing(data,genes_missing)
    missing_ID=[];
    for i in genes_missing
        push!(missing_ID,data[i,:New_ID])
    end 
    return missing_ID
end 


## Calculates the distribution number of randomers that hit each gene at least once. Not very useful. 
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
function gene_coverage_counts(randomers,mRNAfile)
    mRNAdata=CSV.read(mRNAfile,DataFrame);
    mRNAs=read_mRNA_list(mRNAdata)
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
function dnaarray2stringarray(array)
    array2=[]
    for i=1:length(array)
        push!(array2,convert(String,array[i]))
    end
    return array2
end
function stringarray2dnaarray(array)
    array2=[]
    for i=1:length(array)
        push!(array2,LongDNASeq(array[i]))
    end 
    return array2
end 

function pooled_degeneracy(pool)
    deg=0;
    for randomer in pool
        deg +=degeneracy(randomer;uselog2=false);
    end
    return deg
end 
 ########################################################################
 # #CODE FOR CREATING RANDOMER POOLS And Saving them to files
 ########################################################################

 #Quick Test#
 # data=run_NSR_RL(pool_size=3,nsims=50)
 #CSV.write("NSR_RL_S_mutans_4_12_21",data)
function run_NSR_RL(;species="S_mutans",randomerlen=6,pool_size=25,all_bases = dna"AGCTMRWSYKVHDBN",kwargs...)
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
    hitcounts=zeros(numgenes)
    times=[]
    num_blocked_sites=zeros(pool_size)
    println("Initialization Successful")
    for i = 1:pool_size
        if i > 1   
            push!(blocking_sites,LongDNASeq(randomers[i-1]));
            push!(blocking_sites,LongDNASeq(reverse_complement(randomers[i-1]))) ; ## avoid adding randomers that are reverse complements of each other
            blocking_sites=unique(decompress_pool(blocking_sites));
        end  
        num_blocked_sites[i]=length(blocking_sites)
        randomerstats= @timed rollout_rt(all_ns,blocking_sites,mRNAs; hitcounts=hitcounts, simulate=simulate_random_mRNA_uniform_coverage,kwargs...); #create new randomer and add to pool
        push!(randomers,randomerstats.value)
        push!(times,randomerstats.time)
        final_score,hitcounts,dist_score,uniformity_weight,total_hits_weight=mRNA_uniform_coverage(randomers[i],mRNAs; hitcounts=hitcounts, overwrite_hitcounts=true,kwargs...); #officially add hit positions from the new randomer 
        #println("Iteration $i")
    end 
    final_score,hitcounts,dist_score,uniformity_weight,total_hits_weight =mRNA_uniform_coverage(dna"-"^randomerlen,mRNAs; hitcounts=hitcounts, overwrite_hitcounts=false,kwargs...);
    genes_hit=length(findall(x ->x>0,hitcounts))
    genes_missing=findall(x -> x==0,hitcounts)
    total_hits=sum(hitcounts);
    missing_IDs=join(ID_genes_missing(mRNA_DF,genes_missing),",")
    totaldeg=pooled_degeneracy(randomers)
    total_time=sum(times)
    string_time=join(times,",")
    string_num_blocked_sites=join(num_blocked_sites,",")
    string_randomers=join(dnaarray2stringarray(randomers),",")
    data=DataFrame(poolsize=pool_size, primer_length=randomerlen, uniformity_weight=uniformity_weight,hits_weight=total_hits_weight,species=species,num_genes=numgenes,
        transcriptome_size=transcriptome_size,mRNA_GC=mRNA_GC,blocked_GC=rRNA_GC,randomers=string_randomers,degeneracy=totaldeg,num_genes_hit=genes_hit,
        genes_missing=missing_IDs,total_hits=total_hits,coverage_uniformity_score=dist_score,individual_times=string_time,time=total_time,num_blocked_sites=string_num_blocked_sites)
    return data 
end 

function benchmark_NSR_RL()
    DF1=DataFrame()
    nexps=1;
    uniformity_weights=[0]
    total_hits_weights=[100]
    species="S_mutans"
    filename="./S_mutans_NSR_Alt_Valid_Randomer_Test_Timing_6_10_21.csv"
    for i = 1:nexps
        data=run_NSR_RL(;species=species,pool_size=3,nsims=100,uniformity_weight=uniformity_weights[i],total_hits_weight=total_hits_weights[i])
        DF1=vcat(DF1,data);
        CSV.write(filename,DF1)
    end 
end


##Code for analyzing previously made primer pools. stores the same data that run_NSR_RL stores 
function analyze_NSR_Pool(randomers;species="S_mutans",mRNAfile="./SMU_UA159_Genes.csv",rRNAfile="./SMU_UA159_rRNA_tRNA.csv",kwargs...)
    mRNAdata=CSV.read(mRNAfile,DataFrame)
    rRNAdata=CSV.read(rRNAfile,DataFrame)
    mRNAs=read_mRNA_list(mRNAdata);
    rRNAs=read_rRNA_tRNA_list(rRNAdata);
    numgenes=length(mRNAs)
    randomerlen=length(randomers[1]);
    blocking_sites=read_blocking_sites(rRNAdata,randomerlen=randomerlen);
    mRNAbp,mRNA_GC=calculate_transcriptome_stats(mRNAs);
    rRNAbp,rRNA_GC=calculate_transcriptome_stats(rRNAs);
    transcriptome_size=mRNAbp+rRNAbp;
    hitcounts=zeros(numgenes)
    dist_score=0
    for i=1:length(randomers)
        final_score,hitcounts,dist_score,uniformity_weight,total_hits_weight=mRNA_uniform_coverage(randomers[i],mRNAs; hitcounts=hitcounts, overwrite_hitcounts=true,kwargs...); #officially add hit positions from the new randomer 
    end 
    genes_hit=length(filter(x->x>0,hitcounts));
    genes_missing=findall(x->x==0,hitcounts);
    total_hits=sum(hitcounts);
    missing_IDs=join(ID_genes_missing(mRNAdata,genes_missing),",")
    totaldeg=pooled_degeneracy(randomers)
    total_time=0
    string_time="N/A"
    string_randomers=join(dnaarray2stringarray(randomers),",")
    uniformity_weight=0
    hits_weight=0
    pool_size=length(randomers);
    data=DataFrame(poolsize=pool_size, primer_length=randomerlen, uniformity_weight=uniformity_weight,hits_weight=hits_weight,species=species,num_genes=numgenes,
        transcriptome_size=transcriptome_size,mRNA_GC=mRNA_GC,blocked_GC=rRNA_GC,randomers=string_randomers,degeneracy=totaldeg,num_genes_hit=genes_hit,
        genes_missing=missing_IDs,total_hits=total_hits,uniformity_score=dist_score,individual_times=string_time,time=total_time)
    return data
end 

function analyze_all_genomes_NSR_Pools(species_list)
    DF1=DataFrame()
    n=length(species_list)
    for i = 1:n
        randomers=stringarray2dnaarray(CSV.read(string("./Brute Force Pools V2 Compressed/",species_list[i],"_NSR_Brute_Force_Compressed.csv"),DataFrame)[:Sequence])
        mRNAfile=string("./Genome Data/",species_list[i],"_Genes.csv")
        rRNAfile=string("./Genome Data/",species_list[i],"_rRNA_tRNA.csv")
        data=analyze_NSR_Pool(randomers;species=species_list[i],mRNAfile=mRNAfile,rRNAfile=rRNAfile)
        DF1=vcat(DF1,data);
        pct_complete=i/n*100;
        println("$pct_complete % Complete")
        CSV.write("./Brute_Force_Pool_Compressed_Analysis_V2.csv",DF1)
    end 
    
    
end 
function analyze_NSR_Pool_gene_coverage(species_name)
    randomers=stringarray2dnaarray(CSV.read(string("./Brute Force Pools V2/",species_name,"_NSR_Brute_Force_all_randomers.csv"),DataFrame)[:Sequence])
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
    cumulative_hits=zeros(nrandomers)
    cumulative_dist_score=zeros(nrandomers)
    hitcounts=zeros(length(mRNAs))
    randomer_num=1:nrandomers;
    str_randomers=dnaarray2stringarray(randomers);
    cumulative_genes_hit=zeros(nrandomers);
    for i=1:nrandomers
        final_score,hitcounts,dist_score,uniformity_weight,total_hits_weight=mRNA_uniform_coverage(randomers[i],mRNAs; hitcounts=hitcounts, overwrite_hitcounts=true,kwargs...);
        cumulative_hits[i]=sum(hitcounts)
        cumulative_dist_score[i]=dist_score;
        cumulative_genes_hit[i]=sum(hitcounts.>0)
    end 
    data=DataFrame(run=randomer_num,randomer=str_randomers,cumulative_genes_hit=cumulative_genes_hit,cumulative_hits=cumulative_hits,cumulative_dist_score=cumulative_dist_score)
    CSV.write(outputfile,data)
end 

#data=run_NSR_RL(;nsims=100,poolsize=1)
#CSV.write("NSR_RL_Prelims_Genome_Test_4_18_21.csv",data)
#speciesdata=CSV.read("NSR_RL_Prelims_Genome_Test_4_18_21.csv",DataFrame)
#SMUrandomers=stringarray2dnaarray(split(speciesdata[3,10],","));
#countdata=gene_coverage_counts(SMUrandomers,"./E_coli_Genes_4_15_21.csv")
#CSV.write("Ecoli_Coverage_Counts_4_18_21.csv",countdata)

#species_list=CSV.read("./Genome Data/Species_List.csv")[:Species]
#analyze_NSR_Pool_gene_coverage.(species_list)

#species_list=CSV.read("./Genome Data/Species_List.csv")[:Species]
#analyze_all_genomes_NSR_Pools(species_list)
#randomers=stringarray2dnaarray(CSV.read("./Brute Force Pools V2 Compressed/S_mutans_NSR_Brute_Force_Compressed.csv",DataFrame)[:Sequence])
#countdata=gene_coverage_counts(randomers,"./Genome Data/S_mutans_Genes.csv")
#CSV.write("S_mutans_NSR_Comparison_2_gene_coverage_5-19-21.csv",countdata)
#data=CSV.read("S_mutans_Add_Best_Sim_Option_Test_2.csv",DataFrame)
#randomers=stringarray2dnaarray(split(data[2,10],","));
#analyze_cumulative_NSR_Pool(randomers,"S_mutans";outputfile="./S_mutans_Cumulative_NSR_Pool_Best_Sim_Option_6_8_21_Uniformity.csv")
#data=CSV.read("S_mutans_NSR_Num_Blocked_Sites_Test_6_8_21.csv",DataFrame)
#randomers=stringarray2dnaarray(split(data[1,10],","))
#randomers=filter(x -> degeneracy(x;uselog2=false)>0,randomers)
#randomers=decompress_pool(randomers)
#rcomps=reverse_complement.(randomers)
#blocking_sites=read_blocking_sites(CSV.read("./Genome Data/S_mutans_rRNA_tRNA.csv",DataFrame));
#sites=unique(vcat(blocking_sites,rcomps,randomers))
#benchmark_NSR_RL()
