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
- OligoCompressor: compresses a set of non-degenerate oligos into a smaller set of degenerate ones 
- NSR-RL: creates semi-random hexamers that avoid rRNA when preparing RNA-seq libraries
# Instructions to Use OligoRL 
Each software application is a standalone Julia file containing all of the necessary functions to run the application. Requires installation of  [Julia](https://julialang.org/downloads/). Once Julia is installed. Install the dependent Julia Packages. Navigate to package mode by pressing `]` 

```julia 
add Random
add Statistics
add CSV
add Biosequences
add DataFrames
add XLSX
add ArgParse
add RCall
```

# CutFreeRL 
The CutFreeRL application is located in the `CutFreeRL` folder as `CutFreeRL.jl`. CutFree is available as a web-tool at the [Jensen Lab website](http://jensenlab.net/tools/). We have also created a command line interface for the CutFreeRL implementation.


```
$ julia CutFreeRL.jl -sequence NNNNNNNNNNN -restrictionsites GGTCTC,GGCCGG -nsims 1000
```
 

### Required Arguments 
- `--sequence`,`-s`: Starting DNA sequence that should be blocked from containing restriction sites. To generate a set of barcodes with the highest diversity, start with a string of N's the length of your oligo. 

- `--restrictionsites`,`-r`: Sequences to block from the oligo pools. Separate multiple sequences by commas. Do not include spaces.


### Optional Arugments 
- `--nsims`,`-n =1000`: The number of rollout simulations per action. Decreasing `nsims` results in faster runtimes but potentially decreases solution optimality. 


# OligoCompressor
The OligoCompressor application is located in the `OligoCompressor` folder as `OligoCompressor.jl` To use OligoCompressor on your own data, we have created a command line interface.

``` 
$ julia OligoCompressor.jl --targetpool mysequences.txt 
```

### Required Arguments
- `--targetpool`,`-t`: A text file containing a set of non-degenerate oligos to be compressed. An example file called `test_seqs.text` is located in the `OligoCompressor` folder

### Optional Arguments
- `--output`,`-o`: Output .txt file name. The compressed sequences will print if no input is given
- `--bases`,`-b =dna"AGCTMRWSYKVHDBN"`: Nucleotide codes available to be used. Default is all 15.
 
- `--nsims`,`-n =1000`: Number of rollout simulations per action.




# NSR-RL 

The NSR-RL tool reqiures the use of [R](https://www.rstudio.com/products/rstudio/download/) in order to download the genomic data for your organism of interest. We have created the R package [nsrgenomes](https://github.com/bmdavid2/nsrgenomes) to automate this process.  Please refer to [nsrgenomes](https://github.com/bmdavid2/nsrgenomes) to install the package before you use the NSR-RL tool.

The NSR-RL application is located in the `NSR_RL` folder as `NSR_RL.jl`. To use NSR-RL to design your own NSR primers for any organism of interest, we have created a command line interface.


``` 
$ julia  NSR_RL.jl --species E_coli --accession NC_000913.3 --poolsize 100
```
### Arguments
- `--species`,`-s`: The species name of interest. Ex. "E_coli"
- `--poolsize`,`-p`: The number of NSR primers in the final pool

Either
- `--accession`,`-a`: The genome accession number listed on the NCBI database

Or 

- `--genefile`: The file name for the file  containing the mRNA sequences of the organsim of interest made by the nsrgenomes R package
- `--rRNAtRNAfile`: the file name for the file containing the rRNA and tRNA sequences made by the nsrgenomes R package


### Optional Arguments 
- `--outfile`,`-o`: The file name for the NSR-RL output. Provides a more detailed set of statistics.  
- `--length_primers`,`-l =6`: The desired length of each primer. Default is hexamer. 
- `--bases`,`-b =AGCTMRWSYKVHDBN`: The available base codes for the primers
- `--nsims`,`-n =1000`: The number of rollout simulations per action

The following arguments allow uses to change the weights in the reward function: 
    
$$Reward = \beta_{GenesHit}GenesHit + \beta_{TotalHits}TotalHits + \beta_{Inter}Interuniformity + \beta_{Intra}Intrauniformity$$
This allows users to empiracally tune NSR pools. 
- `--genes_hit_weight =1`: Weight given to hitting each gene at least once in the reward function 
- `--total_hits_weight =1`: Weight given to maximizing the total number of hits 
- `--interuniformity_weight =1`: Weight given to placing hits equally across all genes
- `--intrauniformity_weight =1`: Weight given to placing hits equally within each gene 


# Instructions to Replicate Published Data 
We provide pre-written scripts to generate the data shown in each of the figures. Requires installations of [R](https://www.rstudio.com/products/rstudio/download/) ,  [Julia](https://julialang.org/downloads/), and [Gurobi Optimizer](https://www.gurobi.com) (For replicating the original [CutFree](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5762418/) data only). While not required, it is helpful to have a text editing IDE installed such as VScode or atom. 

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
