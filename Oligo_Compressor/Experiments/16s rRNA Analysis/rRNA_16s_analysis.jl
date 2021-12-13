include("/home/bendavid/rRNA_16s_Analysis/OligoCompressor.jl")

using FASTX
#using Plots


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

function hamming_distance_slider(longseq::LongSequence{DNAAlphabet{4}},shortseq::LongSequence{DNAAlphabet{4}})
    sliderlength=length(longseq)-length(shortseq)+1
    dists=zeros(Int64,sliderlength)
    for i=1:sliderlength
        dists[i]=hamming_distance(longseq[i:i+length(shortseq)-1],shortseq)
    end 
    return dists
end 

seqs=String[]
reader= FASTA.Reader(open("/home/bendavid/rRNA_16s_Analysis/rRNA_16s.fasta","r"))

for record in reader
    push!(seqs,sequence(record))
end 
close(reader)


seqs=LongDNASeq.(seqs)
logics=trues(length(seqs))
for i=1:length(seqs)
    if sum(isambiguous.(seqs[i]))>0
        logics[i]=false
    end 
end 
seqs=seqs[logics]


primer=dna"AGAGTTTGATCCTGGCTCAG"

k=collect(0:length(primer))
n=length(seqs)
counts=zeros(length(k))
idxs=zeros(Int64,n,length(k))
for i=1:length(k)
    for j=1:n
        idx=approxsearchindex(seqs[j],primer,i)
        if idx!=0
            counts[i]+=1
            idxs[j,i]=idx
        end 
    end 
end 

pcts=counts./n
#plot(k,pcts)
#histogram(idxs[:,5])
#histogram(idxs[:,6])
#histogram(idxs[:,7])

mismatches=4
valids=idxs[:,mismatches+1].>0 .& idxs[:,mismatches+1].<2
valid_seqs=seqs[valids]


seqs_100=[dna"-" for i=1:n]

for i=1:n
    seqs_100[i]=seqs[i][1:100]
end 


valid_primers=LongDNASeq[]
for i=1:n
    dists=hamming_distance_slider(seqs_100[i],primer)
    val,idx=findmin(dists)
    if val <=4
        push!(valid_primers,seqs_100[i][idx:idx+length(primer)-1])
    end 
end 

valid_primers
unique(valid_primers)

compressed_primers=oligo_pool_compressor(unique(valid_primers))


DF1=DataFrame(sequence=unique(valid_primers))
DF2=DataFrame(sequence=compressed_primers)

CSV.write("./rRNA16s_uncompressed.csv",DF1)
CSV.write("./rRNA16s_compressed.csv",DF2)


