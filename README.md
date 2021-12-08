# <img src="./OligoRL.png" width="30%" align="right" /> OligoRL
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/bmdavid2/OligoRL/blob/main/LICENSE)
# Contents 
[Description](#description) \
[Instructions to Use OligoRL](#instructions-to-use-oligorl) \
[CutFreeRL](#cutfreerl) \
[OligoCompressor](#oligocompressor) \
[NSR-RL](#nsr-rl)\
[Instructions to Replicate Published Data](#instructions-to-replicate-published-data) 

# Description 
A Reinforcemnt Learning Framework for Pooled Oligonucleotide Design. Solve oligo design problems with complex constraints. 

OligoRL has three associated tools: 
- CutFreeRL: finds degenerate barcodes that lack specified restriction enzyme sites
- OligoCompressor: compresses a set of individual oligos into a smaller set of degenerate ones 
- NSR-RL: creates semi-random hexamers that avoid rRNA when preparing RNA-seq libraries
# Instructions to Use OligoRL 
Each software application exists as a standalone Julia file containing all of the necessary functions to run the application. Requires installation of  [Julia](https://julialang.org/downloads/). Once Julia is installed. Install the dependent Julia Packages. Navigate to package mode by pressing `]` 

```julia 
add Random
add Statistics
add CSV
add Biosequences
add DataFrames
add XLSX
add ArgParse
```

# CutFreeRL 
The CutFreeRL application is located in the `CutFreeRL` folder as `CutFreeRL.jl`. Open the file and run the script containing the CutFreeRL functions to use the application.

Users can design degenerate barcodes that lack specified restriction enzyme sites using the `cutfree_rollout` function. `cutfree_rollout` returns a single degenerate randomer.


```julia
    cutfree_rollout(bases, sites; simulate=simulate_random, kwargs...)
```
 

### Arguments 
- `bases`: An array of allowed bases at each postion, usually all 15 degenerate bases. Option is given because some companies restrict which degernate bases are allowed. This is also how to specify the desired randomer length.
#### Example
```julia 
    bases=[dna"AGCTMRWSYKVHDBN" for i =1:10] #Make a 10-mer using all codes
```
- `sites`: An array of restriciton enzyme recognition sequences to be blocked in the random barcode. Encoded as 5' -> 3'.
#### Example
```julia
    sites=[dna"GGTCTC",dna"GAATTC",dna"AAGCTT"] #BsaI,EcoRI,HindIII
```
- `simulate`: The choice of policy for rollout simulations. `simulate_random` will use a random rollout policy. `simulate_greedy` will use a greedy 1 step lookahead policy. 

### Optional Keyword Arugments 
- `nsims=1000`: The number of rollout simulations per action. Decreasing `nsims` results in faster runtimes but potentially decreases solution optimality. 

### Running CutFreeRL with your own data
```julia 
randomer_length=10;
bases=[dna"AGCTMRWSYKVHDBN" for i =1:randomer_length]
sites=[dna"GGTCTC",dna"GAATTC",dna"AAGCTT"] #BsaI,EcoRI,HindIII
randomer=cutfree_rollout(bases,sites; simulate=simulate_random,nsims=1000)


```

# OligoCompressor
The OligoCompressor application is located in the `OligoCompressor` folder as `OligoCompressor.jl` To use OligoCompressor on your own data, we have created a command line interface.

``` 
$ julia OligoCompressor.jl -targetpool mysequences.txt 
```

### Required Arguments
- `targetpool`: A text file containing a set of non-degenerate oligos to be compressed. An example file called `test_seqs.text` is located in the `OligoCompressor` folder

### Optional Arguments
- `nucleotides=dna"AGCTMRWSYKVHDBN"`: Nucleotide codes available to be used. Default is all 15.
### Optional Keyword Arguments 
- `nsims=1000`: Number of rollout simulations per action.

### Running OligoCompressor with your own data
```julia 
    uncompressed_pool=LongDNASeq.(CSV.read("mydata.csv",DataFrame)[:Sequence]) 
    #Sequences stored in a .csv file under a header called "Sequence"
    compressed_pool=oligo_pool_compressor(uncompressed_pool,nucleotides=dna"AGCTMRWSYKVHDBN";nsims=100)
```


# NSR-RL 
The NSR-RL application is located in the `NSR_RL` folder as `NSR_RL.jl`. Open the file and run the script containing the NSR-RL functions to use the application.

Users can design degenerate barcodes that lack specified restriction enzyme sites using the `NSR_RL` function. `NSR_RL` returns data frame containing the designed NSR primers along with other information about the run.

```julia 
    NSR_RL(;species="S_mutans",randomerlen=6,pool_size=25,nucleotides=dna"AGCTMRWSYKVHDBN",kwargs...)
```
### Arguments
- `species="S_mutans"`: The species for which the primers will be designed. Used to extract gene and rRNA_tRNA data from the Genome Data folder in the current directory [see Formatting Genomes for NSR RL](#formatting-genomes-for-nsr-rl)
- `randomerlen=6`: The length of the primers to be designed 
- `pool_size=25`: The number of primers to be made 
- `nucleotides=dna"AGCTMRWSYKVHDBN"`: The base set of available codes to be used for designing primers.

### Optional Keyword Arguments 
- `nsims=1000`: The number of rollout simulations per action

    The following arguments allow uses to change the weights in the reward function: 
    
    $$Reward = \beta_{GenesHit}GenesHit + \beta_{TotalHits}TotalHits + \beta_{Inter}Interuniformity + \beta_{Intra}Intrauniformity$$
    This allows users to empiracally tune NSR pools. 
- `genes_hit_weight=1`: Weight given to hitting each gene at least once in the reward function 
- `total_hits_weight=1`: Weight given to maximizing the total number of hits 
- `interuniformity_weight=1`: Weight given to placing hits equally across all genes
- `intrauniformity_weight=1`: Weight given to placing hits equally within each gene 

### Formatting Genomes for NSR RL
`NSR_RL` requires that the user have a formatted copy of the transcriptome and rRNA/tRNA sequences for the organism of interest. We have provided an R script called `Genome_Info_Mining.R` to automate the genome retrival and formatting process. Note that `Genome_Info_Mining.R` has only been tested for prokaryotic genomes. 

1. To retrieve, format, and save a genome of interest, open `Genome_Info_Mining.R` and navigate to **Line 472**. The next few lines should look like this: 

```r
genome_accession= c("NC_004350.2")#Input list of accession numbers separated by comma
x=lapply(genome_accession,new_df_0_creation) #Output list where each element in the list is a list itself. Element 1 in this list is the rRNA/tRNA seqs and element 2 is the mRNA seqs.
species_name="S_mutans"
suffix_gene="_Genes.csv"
suffix_rRNA_tRNA="_rRNA_tRNA.csv"
```
2. Change the genome accession numbers and species name according to the organsim of interest. If the organsim has multiple accession numbers, separate them by comma in genome accession vector. 
3. Run the entire `Genome_Info_Mining.R` script. Two files will be saved to the working directory `species_name_Genes.csv` and `species_name_rRNA_tRNA.csv`. 
4. These two files can be automatically used by `NSR_RL` by putting them in a folder called `Genome Data` that exists in the working directory of `NSR_RL.jl`. When running `NSR_RL`, the function uses the species argument to extract both genome files from the `Genome Data` folder. The function handles the orientation of the genes, such that the designed NSR primers will theoretically hybridize to the proper strand. 


# Instructions to Replicate Published Data 
We provide pre-written scripts to generate the data shown in each of the figures. Requires installations of [R](https://www.rstudio.com/products/rstudio/download/) ,  [Julia](https://julialang.org/downloads/), and [Gurobi Optimizer](https://www.gurobi.com) (For replicating the original [CutFree](https://jensenlab.shinyapps.io/cutfree/) data only). While not required, it is helpful to have a text editing IDE installed such as VScode or atom. 

### 1. Install required R and Julia Packages 
Install R Packages

```r
install.packages(ggplot2)
install.packages(gurobi)
install.packages(dplyr)
install.packages(stringr)
```
Istall Julia Packages. Navigate to package mode by pressing `]`
```julia 
add Random
add Statistics
add CSV
add Biosequences
add DataFrames
add XLSX
```
### 2. Replicate Figure 1 Data
1. Navigate to `Figures -> Figure1`. Open and run the `rebase.R` file. This file designs the experiments and runs the orignal Cutfree algorithm. Information about the experiments and the CutFree benchmarking data are saved to .csv files in the open directory. 
2. Open and run the `Figure1.jl` file. This file extracts the experimental settings used for the original CutFree algorithm and runs CutFreeRL using these settings. Data for plotting and analysis are saved to .csv files in the open directory
3. Open and run the `plot_Figure1.R` file. This file plots the data as it appears in Figure 1 and performs the statistical analyses reported in the paper. 

### 3. Replicate Figure 2 Data
1. Navigate to `Figures -> Figure2`. Open and run the `Figure2.jl` file. This file contains funcitons that design and run the experiments for OligoCompressor and then saves the data from those experiments to .csv files. Along with the data that appears in the figure, this file also runs the comparision experiment for a published NSR pool described in text. The compressed pool is saved as a .csv file in the open directory.  
2. Open and run the `plot_Figure2.R` file. This file plots the data gathered from the julia script as it appears in Figure 3. 

### 4. Replicate Figure 3 Data
1. Navigate to `Figures -> Figure3`. Open and run the `Figure3.jl` file. This file contains the NSR-RL base funcitons along with funcitons that design and run the experiments described in the text. Data from these experiments are saved to .csv files for visualization and analysis. The Genome Data folder contains the transcriptomes and rRNA/tRNA sequences of 25 prokaryotic species used to benchmark NSR-RL. 

2. Open and run the `plot_Figure3.R` file. This file plots the data gathered from the julia script as it appears in Figure 3. 
