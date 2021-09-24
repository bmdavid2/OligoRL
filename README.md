# <img src="./OligoRL.png" width="30%" align="right" /> OligoRL
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/bmdavid2/OligoRL/blob/main/LICENSE)
## Description 
A Reinforcemnt Learning Framework for Pooled Oligonucleotide Design. Solve oligo design problems with complex constraints. 

OligoRL has three associated tools: 
- CutFreeRL: finds degenerate barcodes that lack specified restriction enzyme sites
- OligoCompressor: compresses a set of individual oligos into a smaller set of degenerate ones 
- NSR-RL: creates semi-random hexamers that avoid rRNA when preparing RNA-seq libraries


## Instructions Replicate Published Data 
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
