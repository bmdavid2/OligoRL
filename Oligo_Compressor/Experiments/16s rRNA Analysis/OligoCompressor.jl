import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences, FASTX, XLSX, CPUTime
import BioSequences: iscompatible

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

# Determine whether or not a candidate sequence belongs to any member of the list of sites. We only want seq that belong to sites 
"""
    isvalid(seq,sites)

Check if the degeneracy of seq is less than or equal to the number of sites it occurs in. 
"""
function isvalid(seq, sites)
    counter=0;
    for site in sites
        if occursin(seq,site) 
            counter+=1; 
        end
    end
    if degeneracy(seq;uselog2=false)<=counter # If the degeneracy of the randomer is greater than the number of sites it hits, then it must hit sites that aren't in the list 
        return true
    else 
    return false 
    end
end

"""
    find_valid_randomer(prefix, bases, sites)

Randomly select bases to add the prefix, ensuring at least one site
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

## Create the list of available actions. Reffered to as either ns or bases in other funcitons 

"""
    define_action_space(sites)

Find the available base codes for an oligo, such that only codes that appear in sites at each postion are included.
"""
function define_action_space(sites)
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
    return bases
end 

""" 
    omit_codes(bases)

Helper function for define_action_space. 
"""
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



#Simulator function 
function simulate_random_hit_score(prefix,bases,sites; nsims=1000, horizon=length(bases), kwargs...)
    hit_scores=zeros(nsims)
    true_horizon = min(horizon, length(bases))
    Threads.@threads for i = 1:nsims   #This is where multithreading is incorporated 
        hit_scores[i]= degeneracy(find_valid_randomer(prefix,bases[1:true_horizon],sites);uselog2=false)
    end
    return mean(hit_scores)
end

### MAIN FUNCTION
"""
    oligo_rollout(bases::Array{LongSequence{DNAAlphabet{4}},1},sites::Array{LongSequence{DNAAlphabet{4}},1}; simulate=simulate_random_hit_score, kwargs...)

Design an oligo that captures as many oligos from sites as possible

# Arguments
- `bases`: An array allowed bases at each postion, usually all 15 degenerate bases. Option is given because some companies restrict which degernate bases are allowed.
- `sites`: An array of non-degenerate oligos to be captured by an oligo. 
- `simulate`: The choice of policy for rollout simulations. simulate_random_hit_score will perform a random rollout with reward for capturing sites. 
"""
function oligo_rollout(bases::Array{LongSequence{DNAAlphabet{4}},1},sites::Array{LongSequence{DNAAlphabet{4}},1}; simulate=simulate_random_hit_score, kwargs...)
    n=length(bases)
    randomer= dna"-" ^ n
    # Find the longest oligo and use this as the horizon. (they should all be equal)
    max_len = map(length, sites) |> maximum
    horizon = max_len
    for i= 1:n
        best_base= DNA_Gap
        best_score=0;
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        bases=define_action_space(sites[findall(j -> iscompatible(j[start:i-1],randomer[start:i-1]),sites)])
        if length(bases[i])==1
            best_base=bases[i][1];
            #println("Shortcut Triggered at position $i")
        else 
            for base in bases[i]
                if ~isvalid(randomer[start:i-1] * base, sites)
                # adding this bases creates sequence not in the oligo pool; skip it
                    continue
                end
                score= simulate(randomer[start:i-1] * base, bases[i+1:end], sites; horizon=horizon,kwargs... )

                # Update use this base if:
                #   1. the mean randomer mRNA coverage is higher than the previous best, OR
                #   2. the randomer mRNA coverage is tied with the previous best, but the new base has
                #      higher individual degeneracy           
                if score > best_score || (score == best_score && degeneracy(base) > degeneracy(best_base))
                    best_score = score
                    best_base = base
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
    end

    return randomer
end

"""
    remove_hit_oligos(randomer,sites)

Remove oligos from sites that are captured by randomer. 
"""
function remove_hit_oligos(randomer,sites)
    n=length(sites);
    hit_indicies=[]
    for i=1:n 
        if occursin(sites[i],randomer)
            push!(hit_indicies,i)
        end
    end
    deleteat!(sites,hit_indicies)
    return sites
end



##Read in the oligo pool. Excel file must be just a list of randomers in column a with no header on sheet 1 
function read_oligo_pool(file)
    xf=XLSX.readxlsx(file);
    x=xf[XLSX.sheetnames(xf)[1]]
    sites=[]
    n=length(x[:]);
    poolarray=x[:];
 for i in 1:n
        seq=(LongDNASeq(string(poolarray[i])))
    push!(sites,seq)
    end
    return sites
end
##Read in the oligo pool as a data frame. The oligos must have a header called "Sequence"
function read_oligo_pool(data)
    sites=[];
    for i=1:nrow(data)
        seq=(LongDNASeq(string(data[:Sequence][i])))
        push!(sites,seq)
    end 
    return sites 
end 
 


#Random pool generator for testing 
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
    
    return LongDNASeq.(randpool)      
end 

"""
    decompress_oligo(comp_oligo)

Expand a degnerate oligo into the non-degnerate oligos it defines. 
"""
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

"""
    decompress_pool(compressed_pool)

Expand a pool of degenerate oligos into a single pool of non-degnerate oligos. 
"""

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
function generate_poissrnd_pool(;randomerlen=6,randpoolsize=4,desiredsize=16)
    randpool=[];
    randomer=dna"-"^randomerlen;
    nucleotides=[dna"ATGC",dna"TRYKMSW",dna"BDHV",dna"N"]
    lambda=(desiredsize/randpoolsize)^(1/randomerlen);
    for i = 1:randomerlen
        deg,=poissrnd_degeneracy(lambda);
        bases=shuffle(nucleotides[deg])
        randomer[i]=bases[1];
    end 
    push!(randpool,randomer)
    while length(randpool)<randpoolsize
        randomer=dna"-"^randomerlen
        for i = 1:randomerlen
            deg, =poissrnd_degeneracy(lambda);
            bases=shuffle(nucleotides[deg])
            randomer[i]=bases[1];
        end
        check=true 
        for j = 1:length(randpool)
            if occursin(randomer,randpool[j])
                check=false
                break 
            end 
        end
        if check 
            push!(randpool,randomer)
        end
        if length(randpool)==length(nucleotides[1])^randomerlen
            break 
        end 
    end
    
    return LongDNASeq.(randpool)     
end 

function poissrnd_degeneracy(lambda)
    pdf=zeros(4);
    cdf=zeros(4);
    for i=1:4
        pdf[i]=poisson_pmf(i,lambda)
    end 
    pdf=pdf./sum(pdf)
    for i=1:4
        cdf[i]=sum(pdf[1:i])
    end 
    r1=rand()
    scale=(cdf.-r1)
    for i=1:4
        if scale[i]>0
            return i , cdf
            continue 
        end 
    end
end 
function poisson_pmf(x,lambda)
    Px=(lambda^x*exp(-lambda))/factorial(x)
end 
function norm_degeneracy(lambda)
    pdfun=zeros(4);
    cdfun=zeros(4);
    for i=1:4
        values=cdf(Normal(lambda, .8),[i,i+1])
        pdfun[i]=values[2]-values[1]
    end 
    pdfun=pdfun./sum(pdfun)
    for i=1:4
        cdfun[i]=sum(pdfun[1:i])
    end 
    r1=rand()
    scale=(cdfun.-r1)
    for i=1:4
        if scale[i]>0
            return i ,cdfun
            continue 
        end 
    end
end 

## Running function ##
"""
    oligo_pool_compressor(uncompressed_pool::Array{LongSequence{DNAAlphabet{4}},1},nucleotides::LongSequence{DNAAlphabet{4}}=dna"AGCTMRWSYKVHDBN";kwargs...)

Compress an uncompresssed_pool of non-degenerate oligos into a smaller compresssed_pool of degenerate oligos. 

# Arguments
- `uncompressed_pool`: A set of non-degenerate oligos to be compressed.
- `nucleotides=dna"AGCTMRWSYKVHDBN"`: Base codes available to be used. 
# Common Keyword Arguments 
- `nsims`: number of rollout simulations per action 
"""
function oligo_pool_compressor(uncompressed_pool::Array{LongSequence{DNAAlphabet{4}},1},nucleotides::LongSequence{DNAAlphabet{4}}=dna"AGCTMRWSYKVHDBN";kwargs...)
    sites=unique(uncompressed_pool);
    original_poolsize=length(sites)
    compressed_pool=[]
    randomerlen=length(sites[1]);
    ns= [nucleotides for _ in 1:randomerlen];
    while length(sites)>0

        randomer=oligo_rollout(ns,sites;kwargs...)
        sites=remove_hit_oligos(randomer,sites);
        push!(compressed_pool,randomer);
        #checkpoint=length(compressed_pool);
        #println("Compressed Pool Size: $checkpoint")
    end 
    n=length(compressed_pool)
    failed_sites=[]
    for i = 1:length(compressed_pool)
        if degeneracy(compressed_pool[i];uselog2=false)==0
            push!(failed_sites,i)
        end
    end
    deleteat!(compressed_pool,failed_sites);
    new_poolsize=length(compressed_pool)
    println("Compressed Pool contains $new_poolsize oligos. Original Pool contained $original_poolsize oligos.")
    return compressed_pool
end




function benchmark_oligo_compressor(;all_bases=dna"AGCTMRWSYKVHDBN",kwargs...) 
    nreps=10;
    testsizes=[100];
    poolsize=repeat(testsizes,nreps)
    nruns=length(poolsize)
    run=1:nruns
    time=zeros(nruns)
    compressedsize=zeros(nruns)
    pools=[]
    for i = 1:nruns
        randpool=generate_random_pool(;randomerlen=6, randpoolsize=poolsize[i])
        strpool=join(convert.(String,randpool),",");
        push!(pools,strpool)
        stats= @timed oligo_pool_compressor(randpool,all_bases;kwargs...)
        compressedsize[i]=length(stats.value)
        time[i]=stats.time
        filename="oligo_compressor_4_15_21.csv"
        data=DataFrame(poolsize=poolsize,compressedsize=compressedsize,time=time)
        #CSV.write(filename,data)
    end 
    data=DataFrame(pool=pools,poolsize=poolsize,compressedsize=compressedsize,time=time)
    return data
end
function run_recompression_experiment(;nruns=50,all_bases = dna"AGCTMRWSYKVHDBN",kwargs...)
    sizes=repeat([1,5,10,15,20],10)
    runs=1:nruns
    pools=[]
    poolsize=zeros(nruns)
    decompressedsize=zeros(nruns)
    recompressedsize=zeros(nruns)
    for i = 1:nruns 
        randpool=generate_poissrnd_pool(;randpoolsize=10, desiredsize=sizes[i])
        poolsize[i]=length(randpool)
        strpool=join(convert.(String,randpool),",");
        push!(pools,strpool)
        uncomrandpool=decompress_pool(randpool);
        decompressedsize[i]=length(uncomrandpool);
        recomprandpool=oligo_pool_compressor(uncomrandpool,all_bases;kwargs...)
        recompressedsize[i]=length(recomprandpool);
    end 
    data=DataFrame(A=runs,B=pools,C=poolsize,D=decompressedsize,E=recompressedsize)
    rename!(data,[:runs,:pools,:poolsize,:decompressedsize,:recompressedsize])
    CSV.write("recompression_experiment_3-19-21.csv",data)
end 
#benchmark_oligo_compressor()
function rerun_recompression_experiment(filename;nruns=50,all_bases = dna"AGCTMRWSYKVHDBN",kwargs...)
    sizes=repeat([1,5,10,15,20],10)
    runs=1:nruns

    data=CSV.read("recompression_experiment_3-19-21.csv",DataFrame)
    data[:recompressed100]=0;
    for i = 1:nruns 
        randpool=LongDNASeq.(split(data[i,2],","))

        uncomrandpool=decompress_pool(randpool);
        recomprandpool=oligo_pool_compressor(uncomrandpool,all_bases;kwargs...)
        data[i,:recompressed100]=length(recomprandpool);
    end 
    CSV.write(filename,data)
end 


#rerun_recompression_experiment("recompression_3_19_21_add_nsims100.csv";nsims=100)
#data=benchmark_oligo_compressor(;nsims=100)

function run_compression_experiment()
    orig_pool=LongDNASeq.(CSV.read("./Oligo_Compressor/Experiments/Nature_2018_NSR_Primers.csv",DataFrame)[:Sequence])
    all_bases = dna"AGCTMRWSYKVHDBN"
    new_pool=oligo_pool_compressor(orig_pool,all_bases;nsims=100)
    data=DataFrame(randomers=new_pool)
    orig_len=length(orig_pool)
    new_len=length(new_pool)
    println("Original Length: $orig_len, New Length: $new_len")
    filename="./Oligo_Compressor/Experiments/Nature_2018_NSR_Primers_Compressed.csv"
    CSV.write(filename,data)
end
