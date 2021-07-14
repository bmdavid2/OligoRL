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
        candidates=[];
        array=1:length(bases[i])
        candorder=shuffle(array)
        for j in 1:length(array)
        push!(candidates,bases[i][candorder[j]])
        end
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
function hit_stats(indicies,mRNAs,randomerlen)
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
function randomer_coverage_counts(randomers,mRNAfile)
    mRNAdata=CSV.read(mRNAfile,DataFrame);
    mRNAs=read_mRNA_list(mRNAdata)
    randomer_hits=zeros(length(randomers))
    randomerlen=length(randomers[1]);
    for i in 1:length(randomers)
        counter=0;
        for j = 1:length(mRNAs)
            for k in 1:length(mRNAs[j])-randomerlen
                if occursin(randomers[i],mRNAs[j][k:k+randomerlen-1])
                    counter+=1
                end
            end
        end
        randomer_hits[i]=counter
    end
    data=DataFrame(randomers=randomers,counts=randomer_hits)
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


function NSR_Brute_Force(species_name)
    rRNAfile=string("./Genome Data/",species_name,"_rRNA_tRNA.csv")
    mRNAfile=string("./Genome Data/",species_name,"_Genes.csv")
    sites=read_blocking_sites(CSV.read(rRNAfile,DataFrame));
    randomers=generate_random_pool(;randpoolsize=4^length(sites[1]))
    nsites=length(sites);
    del_indicies=[]
    for i=1:length(randomers)
        for j=1:nsites
            if randomers[i]==sites[j]
                push!(del_indicies,i)
                break 
            end 
        end 
    end 
  
    deleteat!(randomers, del_indicies)
    checkpoint=length(randomers)
    first_del=[];
    for j=1:length(randomers)
        if is_palindrome(randomers[j]) 
            push!(first_del,j)
        end 
    end 
    deleteat!(randomers,first_del)
    checkpoint=length(randomers)
    sec_del=[];
    rc_cands=[]
    rcs=reverse_complement.(randomers)
    for randomer in randomers
        if !isvalid(randomer,rcs)
            cand1=randomer
            cand2=reverse_complement(randomer)
            push!(rc_cands,cand1)
            push!(rc_cands,cand2)
        end
    end 
    rc_cands=unique(rc_cands)
    rcpairs=length(rc_cands)/2;
    println("$rcpairs")
    hit_stats=randomer_coverage_counts(rc_cands,mRNAfile)
    for i=1:2:length(rc_cands)
        v,idx=findmin(hit_stats[i:i+1,2])
        del_cand=rc_cands[i+idx-1]
        push!(sec_del,findall(x->x==del_cand,randomers)[1])
    end 
    sec_del=sort(unique(sec_del))
    deleteat!(randomers,sec_del)
    data=DataFrame(Sequence=randomers);
    outputfile=string("./Brute Force Pools V2/",species_name,"_NSR_Brute_Force_all_randomers.csv")
    CSV.write(outputfile,data)

end 

function record_pool_size(species)
    n=length(species)
    poolsize=zeros(n)
    for i=1:n
        poolfile=string("./Brute Force Pools V2/",species[i],"_NSR_Brute_Force_all_randomers.csv");
        poolsize[i]=nrow(CSV.read(poolfile,DataFrame));
    end 
    DF=DataFrame(species=species,poolsize=poolsize)
    outputfile="./Brute Force Pool Sizes V2.csv"
    CSV.write(outputfile,DF)
end 
species=CSV.read("./Genome Data/Species_List.csv",DataFrame)
NSR_Brute_Force.(species[:Species])
record_pool_size(species[:Species])

