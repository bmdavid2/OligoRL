#################################################################
##########Genome Information Mining, mRNA and rRNA###############
#################################################################

 #install.packages("BiocManager")
 #BiocManager::install("genbankr")
 #install.packages("devtools")
# devtools::install_github("jensenlab/primer3")
#devtools::install_github("mhahsler/rBLAST")

library(genbankr)
library(primer3)
library(kmer)
library(Biostrings)
library(insect)
library(tictoc)
library(plyr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(snow)
library(corrplot)
library(reshape2)
library(phylotools)
library(rentrez)

transcript_sequence_finder=function(genome,gene){ #The purpose of this function is to determine the sequence 5'->3' of the coding strand for a gene (what mRNA will look like)
  seq=as.character(genome@sequence[1])
  seq=unname(seq)
  range=ranges(gene)
  start=as.integer(start(range))
  end=as.integer(end(range))
  if(as.character(gene@strand)== "+"){ #If the gene is on the positive strand then the provided start index is the start of the transcript
    trans_seq=as.character(substr(seq,start,end)) #Add buffer zone to account
  }else if(as.character(gene@strand)== "-"){#If the gene is on the negative strand then the start index is the end of the transcript
    trans_seq=as.character(rc(substr(seq,start,end)))
  }
  return(trans_seq)
}

transcript_plus_buffer_sequence_finder=function(genome,gene,lig_junction=12,total_hybrid=30){ #The purpose of this function is to determine the sequence 5'->3' of the coding strand for a gene (what mRNA will look like)
  seq=as.character(genome@sequence[1])
  seq=unname(seq)
  range=ranges(gene)
  start=as.integer(start(range))
  end=as.integer(end(range))
  flanker=(total_hybrid-lig_junction)/2
  
  if(as.character(gene@strand)== "+"){#If the gene is on the positive strand then the provided start index is the start of the transcript
    if(gene$locus_tag=="SMU_RS09895"){ #Special case where the buffer zone wraps around to the start of the index
      trans_seq1=as.character(substr(seq,start = (start-flanker),stop = (end+8)))
      trans_seq2=as.character(substr(seq,start = 1,stop = 1))# Need to change this to make it general depending on flanker
      trans_seq=paste0(trans_seq1,trans_seq2)
    }
    if(gene$locus_tag=="DLJ52_00005"){ #Special case where the buffer zone wraps around to the start of the index
      trans_seq1=as.character(substr(a,start=2144746-(flanker-1),stop=2144746))
      trans_seq2=as.character(substr(seq,start=start,stop=(end+flanker)))# Need to change this to make it general depending on flanker
      trans_seq=paste0(trans_seq1,trans_seq2)
    }else{
      trans_seq=as.character(substr(seq,start = (start-flanker),stop=(end+flanker))) #Add buffer zone to account
    }
  }else if(as.character(gene@strand)== "-"){#If the gene is on the negative strand then the start index is the end of the transcript
    if(gene$locus_tag=="BTH_I0001"){ #Special case where the buffer zone wraps around to the start of the index
      trans_seq1=as.character(rc(substr(seq,start,stop=(end+flanker))))
      trans_seq2=as.character(rc(substr(seq,start=(3809201-(flanker-1)),stop = 3809201)))
      trans_seq=paste0(trans_seq1,trans_seq2)
    }else{
      trans_seq=as.character(rc(substr(seq,start=(start-flanker),stop=(end+flanker))))
    }
  }
  return(trans_seq)
}


hasher=function(x, n){ #This function returns all subset n-mers from a character string,x
  return(substring(x, 1:(nchar(x) - n + 1), n:nchar(x)))
} 


mutate_sequence=function(string, num = 1, nucleotides = c("A","T","C","G")) { ##Not run but a useful function to have, give it a DNA string and it'll generate all possible mutants with a given number of substitutions
  l_str=str_length(string)
  
  choices <- cross(list(
    cols = combn(seq_len(l_str), num, simplify = F),
    muts = cross(rerun(num, nucleotides)) %>% map(unlist)
  ))
  
  choice_matrix <- 
    map_dfr(choices, as_tibble, .id = "rows") %>% 
    mutate(rows = as.numeric(rows))
  
  seq_matrix <- str_split(rep(string, max(choice_matrix$rows)), "", simplify = T)
  
  seq_matrix[as.matrix(choice_matrix[,1:2])]
  apply(seq_matrix, 1, paste, collapse = "")
}

kmer_freq_calc=function(input_vector,kmer_length=c(4,6),freq_table1,freq_table2){ #MAKE SURE THE FREQ TABLES PASSED MATCH THE kmers you want.This function accepts a vector of ligation junction length k-mers (All targeting the same transcript) and gives frequency scores for each subset kmer of specified length (k1,k2), with the default being 4 and 6
  helper=function(mer,df){ #This function takes a k1-mer sequence, searches the genome's k1-mer frequency table for that particular k1-mer, and returns the number of times that k1-mer appeared in the genome 
    df$x=as.character(df$x)
    index=which(df$x == mer)
    output=df$freq[index]
    return(output)
  }
  df1=freq_table1
  df2=freq_table2
  output=list()
  for(i in 1:length(input_vector)){
    kmer1_hashes=hasher(input_vector[i],kmer_length[1]) #Generates all k1 kmers in a given ligation junction or section
    kmer1_scores=vapply(kmer1_hashes,helper,df=df1,FUN.VALUE = 0) #Finds the frequency score for each k1 kmer in a given ligation junction
    kmer1_total=sum(kmer1_scores) #Sums all frequencies
    kmer1_max=max(kmer1_scores) # Gives frequency of the most common kmer1 in input
    
    kmer2_hashes=hasher(input_vector[i],kmer_length[2]) #Generates all k2 kmers in a given ligation junction
    kmer2_scores=vapply(kmer2_hashes,helper,df=df2,FUN.VALUE = 0) #Finds the frequency score for each k2 kmer in a given ligation junction
    kmer2_total=sum(kmer2_scores) #Sums all frequencies
    kmer2_max=max(kmer2_scores) #Gives frequency of the most common kmer2 in input
    
    bundle=c(kmer1_total,kmer1_max,kmer2_total,kmer2_max) #Combines k1 and k2 scores
    names(bundle)=c("kmer1 Freq. Sum","kmer1 Max","kmer2 Freq. Sum","kmer2 Max") #Assigns names
    output[[i]]=bundle #Stores ligation junction scores in a list.
    
  }
  
  return(output) # Returns an output list where each element is a vector carrying the k1 and k2 kmer frequency scores
  
}

flank_add=function(ligation_junction,rc_mRNA,flank_length){ #This function accepts a ligation junction for a given gene, the rc of the mRNA transcript, and the length of the flank sequence to add
  reg_ex_seq=paste0("[A-Z]{",flank_length,"}",ligation_junction,"[A-Z]{",flank_length,"}")
  output=unlist(regmatches(rc_mRNA,gregexpr(reg_ex_seq,rc_mRNA,perl=TRUE)))[1]
  return(output)
}


homopolymer_filter=function(input_df,input_vector,max_length=4){#Returns the indexes of homopolymer tract violators
  Ahomo=as.character(strrep("A",max_length+1))
  Thomo=as.character(strrep("T",max_length+1))
  Chomo=as.character(strrep("C",max_length+1))
  Ghomo=as.character(strrep("G",max_length+1))
  df=input_df
  sequence_vector=input_vector
  a=grep(Ahomo,sequence_vector,fixed = TRUE)
  t=grep(Thomo,sequence_vector,fixed = TRUE)
  c=grep(Chomo,sequence_vector,fixed = TRUE)
  g=grep(Ghomo,sequence_vector,fixed = TRUE)
  index=unique(c(a,t,c,g))
  df=df[-index,]
  return(df)
}

add_ligation_pair=function(input_df){
  helper1=function(seq,size_junc){
    temp=unlist(strsplit(seq,""))
    donor=temp[((size_junc/2)+1)]
    acceptor=temp[(size_junc/2)] 
    output=paste0(acceptor,donor)
    return(output)
  }
  df=input_df
  lig_junction_seqs=df$Probe_Ligation_Junction_Seq
  size_junc=nchar(lig_junction_seqs[1])
  lig_pairs=vapply(lig_junction_seqs,helper1,size_junc=size_junc,FUN.VALUE = "A")
  df$Ligation_Junction_Pair=unname(lig_pairs)
  return(df)
}

add_ligation_pair_efficiency=function(input_df,ligase_data){
  helper2=function(lig_pair,enzyme_data){
    df=enzyme_data
    df$Lig_Pair=as.character(df$Lig_Pair)
    index=which(df$Lig_Pair == lig_pair)
    output1=df$Relative_Yeild[index]
    return(output1)
  }
  df=input_df
  ligation_pairs=df$Ligation_Junction_Pair
  lig_pair_eff=vapply(ligation_pairs,helper2,enzyme_data=ligase_data,FUN.VALUE = 1)
  df$Ligation_Pair_Efficiency=(unname(lig_pair_eff))
  return(df)
}

add_ligation_pair_rxn_vel=function(input_df,ligase_data){
  helper3=function(lig_pair,enzyme_data){
    df=enzyme_data
    df$Lig_Pair=as.character(df$Lig_Pair)
    index=which(df$Lig_Pair == lig_pair)
    output1=df$Vo_Eo[index]
    return(output1)
  }
  df=input_df
  ligation_pairs=df$Ligation_Junction_Pair
  lig_pair_vel=vapply(ligation_pairs,helper3,enzyme_data=ligase_data,FUN.VALUE = 1)
  df$Ligation_Pair_Vo_Eo=unname(lig_pair_vel)
  return(df)
}

add_probe_binding_location=function(input_df){
  transcript_location=function(hybridization_sequence,mRNA){#Input the probe ligation junction and the target mRNA transcript sequence
    lig_rc=rc(hybridization_sequence)
    start_index=regexpr(lig_rc,mRNA)[1]
    lig_index=start_index+((nchar(lig_rc)/2)-1)
    trans_length=nchar(mRNA)
    lig_location=lig_index/trans_length
    return(lig_location)
  }
  df=input_df
  ProbeBinding_Location=c()
  for(i in 1:length(df$New_ID)){
    ProbeBinding_Location[i]=transcript_location(df$Probe_Ligation_Junction_Seq[i],df$mRNA_seq[i])[1] #Can probably do this with vapply to speed up
  }
  
  df$Probe_Binding_Location=100*ProbeBinding_Location
  return(df)
}

probe_end_Tm_diff=function(input_df){
  df=input_df
  df$Tm_ends_diff=abs(df$Hybrid_3prime_Tm-df$Hybrid_5prime_Tm)
  return(df)
}

import_genome=function(gb_accession){
  id=GBAccession(gb_accession) #UA159 S. mutans Genome
  genome=readGenBank(id)
  return(genome)
}
initial_df=function(genome,sample_size,sample=TRUE){
  genes=genome@genes
  ORFs=genome@genes #Remeber to go back and target unannotated genes like ComS
  ncRNA=genome@other_features
  rRNA=ncRNA[ncRNA$type=="rRNA"]
  tRNA=ncRNA[ncRNA$type=="tRNA"]
  sneaky_t=ORFs$locus_tag%in% tRNA$locus_tag #need to eliminate tRNAs which are also included in ORFS
  ORFs=ORFs[!sneaky_t]
  sneaky_r=ORFs$locus_tag %in% rRNA$locus_tag #need to eliminate rRNAs which are also included in ORFS
  ORFs=ORFs[!sneaky_r]
  if(sample == TRUE){
    rand_index=sample(1:length(ORFs),sample_size,replace=FALSE)
    samp_genes=ORFs[rand_index]
  }else{
    samp_genes=ORFs
  }
  species_info=as.character(samp_genes@seqnames)
  new_gene_info=samp_genes$locus_tag
  starts=start(ranges(samp_genes))
  ends=end(ranges(samp_genes))
  strand=as.character(samp_genes@strand)
  samp_trans=c()
  samp_buffer_trans=c()
  for(i in 1:length(samp_genes)){
    samp_trans[i]=transcript_sequence_finder(genome,samp_genes[i])
    samp_buffer_trans[i]=transcript_plus_buffer_sequence_finder(genome,samp_genes[i],lig_junction = 12,total_hybrid = 30)
  }
  df=data.frame(Species=species_info,New_ID=new_gene_info,Start=starts,End=ends,Strand=strand,mRNA_seq=samp_trans,Buffered_mRNA_seq=samp_buffer_trans,stringsAsFactors = FALSE)
  df$New_ID=as.factor(df$New_ID)
  df$rc_mRNA_seq=rc(df$mRNA_seq)
  df$rc_Buffered_mRNA_Seq=rc(df$Buffered_mRNA_seq)
  df$Length=df$End-df$Start
  return(df)
}

add_probe_seqs=function(input_df,lig_junc_length=12,total_hybridization_length=30){
  library(snow)
  hasher=function(x, n){ #This function returns all subset n-mers from a character string,x
    return(substring(x, 1:(nchar(x) - n + 1), n:nchar(x)))
  } 
  df=input_df
  all_rc_trans=df$rc_mRNA_seq
  all_rc_buffer_trans=df$rc_Buffered_mRNA_Seq
  
  num_threads=detectCores()-1
  cl=makeCluster(num_threads)
  
  target_seqs = parLapply(cl, as.list(all_rc_buffer_trans), hasher, n=total_hybridization_length)
  num_replications=unlist(parLapply(cl,target_seqs,length))
  stopCluster(cl)
  df$reps=num_replications
  df2=df[rep(row.names(df), df$reps), 1:(dim(df)[2]-1)]
  df2$Probe_Hybrid_Seq=as.character(unlist(target_seqs))
  
  unique_indices=!(duplicated(df2$Probe_Hybrid_Seq)|duplicated(df2$Probe_Hybrid_Seq,fromLast=TRUE))
  df2=df2[unique_indices,] ##Eliminate non-unique target sites
  
  hybrid_region=total_hybridization_length
  lig_region=lig_junc_length
  flanker=(hybrid_region-lig_region)/2
  lig_junc_seqs=c()
  first_halfs=c()
  last_halfs=c()
  
  for(i in 1:length(df2$Probe_Hybrid_Seq)){
    lig_junc_seqs[i]=substr(df2$Probe_Hybrid_Seq[i],flanker+1,hybrid_region-flanker)
    first_halfs[i]=substr(df2$Probe_Hybrid_Seq[i],1,(hybrid_region/2))
    last_halfs[i]=substr(df2$Probe_Hybrid_Seq[i],(hybrid_region/2)+1,hybrid_region)
  }
  
  df2$Probe_Hybrid3prime_Seq=first_halfs
  df2$Probe_Hybrid5prime_Seq=last_halfs
  df2$Probe_Ligation_Junction_Seq=lig_junc_seqs
  return(df2)
  
}


unique_filter=function(input_df,input_vector){
  df=input_df
  unique_indices=!(duplicated(input_vector)|duplicated(input_vector,fromLast=TRUE))
  df=df[unique_indices,] ##Eliminate non-unique target sites
  return(df)
}


add_Tm=function(input_df){
  df=input_df
  df$Probe_Ligation_Junction_Tm=primer3::calculate_tm(df$Probe_Ligation_Junction_Seq)
  df$Probe_Hybrid3prime_Tm=primer3::calculate_tm(df$Probe_Hybrid3prime_Seq)
  df$Probe_Hybrid5prime_Tm=primer3::calculate_tm(df$Probe_Hybrid5prime_Seq)
  df$Probe_Hybrid_Tm=primer3::calculate_tm(df$Probe_Hybrid_Seq)
  df$Probe_End_Tm_Diff=abs(df$Probe_Hybrid3prime_Tm-df$Probe_Hybrid5prime_Tm)
  
  return(df)
}

Tm_diff_filter=function(input_df,max_Tm_diff=5){
  df=input_df
  df=df[df$Probe_End_Tm_Diff <= max_Tm_diff,]
  return(df)
  
}

Tm_hybridization_filter=function(input_df,min_probe_Tm=25){
  df=input_df
  df=df[df$Probe_Hybrid3prime_Tm >= min_probe_Tm,]
  df=df[df$Probe_Hybrid5prime_Tm >= min_probe_Tm,]
  return(df)
}

hairpin_filter=function(input_df){
  df=input_df
  df$Probe_Hybrid_3prime_Hairpin_Max_Temp=primer3::calculate_hairpin(df$Probe_Hybrid3prime_Seq)$temp
  df$Probe_Hybrid_5prime_Hairpin_Max_Temp=primer3::calculate_hairpin(df$Probe_Hybrid5prime_Seq)$temp
  df=df[df$Probe_Hybrid_3prime_Hairpin_Max_Temp == 0,]
  df=df[df$Probe_Hybrid_5prime_Hairpin_Max_Temp == 0,]
  return(df)
}

rRNA_tRNA_kmer_database=function(genome,kmer_length=30){
  ncRNA=genome@other_features
  rRNA=ncRNA[ncRNA$type=="rRNA"]
  tRNA=ncRNA[ncRNA$type=="tRNA"]
  rRNA_transcripts=c()
  tRNA_transcripts=c()
  for(i in 1:length(rRNA)){
    rRNA_transcripts[i]=transcript_sequence_finder(genome,rRNA[i])
  }
  
  for(i in 1:length(tRNA)){
    tRNA_transcripts[i]=transcript_sequence_finder(genome,tRNA[i])
  }
  rRNA_transcripts=unique(rRNA_transcripts)
  tRNA_transcripts=unique(tRNA_transcripts)
  rc_rRNA_transcripts=rc(rRNA_transcripts) 
  rc_tRNA_transcripts=rc(tRNA_transcripts)
  hashed_rc_rRNA=list()
  hashed_rc_tRNA=list()
  for(i in 1:length(rc_rRNA_transcripts)){
    hashed_rc_rRNA[[i]]=hasher(rc_rRNA_transcripts[i],kmer_length)
  }
  for(i in 1:length(rc_tRNA_transcripts)){
    hashed_rc_tRNA[[i]]=hasher(rc_tRNA_transcripts[i],kmer_length)
  }
  no_go=unique(c(unlist(hashed_rc_rRNA),unlist(hashed_rc_tRNA))) #These are probe hybridization region sequences that need to be avoided in order to avoid capturing rRNA and tRNA
}


rRNA_tRNA_filter=function(input_df,maximum,database){
  df=input_df
  max_allowed_mismatches=maximum
  probe_hybridization_length=nchar(df$Probe_Hybrid3prime_Seq)[1]
  num_threads=detectCores()-1
  cl=makeCluster(num_threads)
  fuzzy_filter_results1= unlist(parLapply(cl, as.list(df$Probe_Hybrid3prime_Seq), fuzzy_filter,database=database,maximum=max_allowed_mismatches)) #Applies the fuzzy filter to all sampled rc transcript sequences to make sure they avoid rRNA and tRNA, for 50 genes and paralellization this takes ~385 seconds
  fuzzy_filter_results2= unlist(parLapply(cl, as.list(df$Probe_Hybrid5prime_Seq), fuzzy_filter,database=database,maximum=max_allowed_mismatches))
  stopCluster(cl)
  fuzzy_filter_results=fuzzy_filter_results1&fuzzy_filter_results2
  df=df[fuzzy_filter_results,]
  return(df)
}

fuzzy_filter=function(query,database,maximum){ #This function takes a query_set (array of candidate ligation junction sequences) and returns a logical array which indicates whether the element passed the filter or not, provided a set number of mismatches to map from the query to the database
  x=agrep(query, database, max.distance=list(all=maximum, insertions=0, deletions=0, substitutions=maximum)) # Find different function that stops after detecting a single match, look at limiting output
  if(length(x)==0){ #agrep returns members of the database which can be matched to the query string. So if nothing is returned, the query passes
    output=TRUE
  }else{
    output=FALSE
  }
  return(output)
}

generate_kmer_freq_tables=function(kmers=c(4,6),genome=genome){
  k1=kmers[1]
  k2=kmers[2]
  genome_seq=as.character(genome@sequence[[1]])
  k1mers=hasher(genome_seq,k1)
  k2mers=hasher(genome_seq,k2)
  k1mer_df=plyr::count(k1mers)
  total_k1mer=length(k1mers)
  k2mer_df=plyr::count(k2mers)
  total_k2mer=length(k2mers)
  
  k1mer_plot<-ggplot(data=k1mer_df, aes(x=x, y=freq)) +
    geom_bar(stat="identity")+labs(title="k-mer Spectra for UA159 (k=4)", y="Frequency", x="k-mer")+theme(text = element_text(size=5.5),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.25))
  k2mer_plot<-ggplot(data=k2mer_df, aes(x=x,y=freq)) +
    geom_bar(stat="identity")+labs(title="k-mer Spectra for UA159 (k=6)", y="Frequency", x="k-mer")+theme(text = element_text(size=5.5),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.25))
  
  output=list()
  output[[1]]=k1mer_df
  output[[2]]=k1mer_plot
  output[[3]]=k2mer_df
  output[[4]]=k2mer_plot
  return(output)
  
}

new_df_0_creation=function(accession_number){
  
  helper=function(list,element_index){
    return(list[[element_index]])
  }
  
    genome_accession=accession_number 
    genome=import_genome(genome_accession)
    ncRNA=genome@other_features
    rRNA=ncRNA[ncRNA$type=="rRNA"]
    tRNA=ncRNA[ncRNA$type=="tRNA"]
    rRNA_transcripts=c()
    tRNA_transcripts=c()
    if (length(rRNA)>0){
      for(j in 1:length(rRNA)){
        rRNA_transcripts[j]=transcript_sequence_finder(genome,rRNA[j])
      }
    }
    if (length(tRNA)>0){
      for(k in 1:length(tRNA)){
        tRNA_transcripts[k]=transcript_sequence_finder(genome,tRNA[k])
      }
    }
    
    rRNA_transcripts=unique(rRNA_transcripts)
    tRNA_transcripts=unique(tRNA_transcripts)
    bad_guys=c(rRNA_transcripts,tRNA_transcripts)
    output1=data.frame(Sequence=bad_guys)
    df_0=initial_df(genome = genome,sample = FALSE)
    output2=df_0
    output=list(output1,output2)
    return(output)
  
}

save_single_chromosome=function(genome_accession,species_name){
  x=lapply(genome_accession,new_df_0_creation) #Output list where each element in the list is a list itself. Element 1 in this list is the rRNA/tRNA seqs and element 2 is the mRNA seqs.
  suffix_gene="_Genes.csv"
  suffix_rRNA_tRNA="_rRNA_tRNA.csv"
  genefile=paste(species_name,suffix_gene,sep="")
  rRNA_tRNAfile=paste(species_name,suffix_rRNA_tRNA,sep="")
  write.csv(x[[1]][[1]],rRNA_tRNAfile)
  write.csv(x[[1]][[2]],genefile)
}

#########################
genome_accession= c("NC_004350.2")#Input list of accession numbers separated by comma
x=lapply(genome_accession,new_df_0_creation) #Output list where each element in the list is a list itself. Element 1 in this list is the rRNA/tRNA seqs and element 2 is the mRNA seqs.
species_name="S_mutans"
suffix_gene="_Genes.csv"
suffix_rRNA_tRNA="_rRNA_tRNA.csv"
genefile=paste(species_name,suffix_gene,sep="")
rRNA_tRNAfile=paste(species_name,suffix_rRNA_tRNA,sep="")
genedata=data.frame()
rRNA_tRNAdata=data.frame()
for (i in 1:length(genome_accession)){
  genedata =rbind(genedata,x[[i]][[2]])
  rRNA_tRNAdata=rbind(rRNA_tRNAdata,x[[i]][[1]])
}
write.csv(rRNA_tRNAdata,rRNA_tRNAfile)
write.csv(genedata,genefile)

















