setwd("~/Documents/Jensen Lab/OligoRL/NSR Data Analysis")
library(ggplot2)
library(stringr)
bfdata <- read.csv("./Brute_Force_Pool_Compressed_Analysis_V2.csv")
names(bfdata)
bfdata$Diff_GC <- bfdata$mRNA_GC-bfdata$blocked_GC
bfdata$hits_per_gene <- bfdata$total_hits/bfdata$num_genes
ggplot(bfdata,aes(x=Diff_GC,y=uniformity_score))+
          geom_point()+
          xlab("GC% Difference Between mRNA and rRNA/tRNA")+
          ylab("Pool Uniformity Score")+
  theme_classic()
ggplot(bfdata,aes(x=Diff_GC,y=hits_per_gene))+
  geom_point()+
  xlab("GC% Difference Between mRNA and rRNA/tRNA")+
  ylab("Average Hits/Gene")+
  theme_classic()
ggplot(bfdata,aes(x=poolsize,y=hits_per_gene))+
  geom_point()+
  xlab("Pool Size")+
  ylab("Average Hits/Gene")+
  theme_classic()
ggplot(bfdata,aes(x=transcriptome_size,y=mRNA_GC))+
  geom_point()+
  xlab("Exome Size")+
  ylab("mRNA GC%")+
  theme_classic()
ggplot(bfdata,aes(x=transcriptome_size,y=hits_per_gene))+
  geom_point()+
  xlab("Exome Size")+
  ylab("Average Hits/Gene")+
  theme_classic()
####################################################3


SMU_cumulative <- read.csv("S_mutans_Cumulative_NSR_Pool_Data_5_18_21_3.csv")
SMU_cumulative$failed <- grepl("-",SMU_cumulative$randomer)
SMU_cumulative <- SMU_cumulative[which(SMU_cumulative$failed==FALSE), ]
new_run=1:nrow(SMU_cumulative)
SMU_cumulative$new_run <- new_run
PlotA <- ggplot(SMU_cumulative,aes(x=new_run,y=cumulative_dist_score))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=0.337,colour="black")+
  xlab("Pool Size")+
  ylab("DUlq Score")+
  theme_classic()
PlotB <- ggplot(SMU_cumulative,aes(x=new_run,y=cumulative_hits))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=412260,colour="black")+
  xlab("Pool Size")+
  ylab("Total Hits")+
  theme_classic()
PlotA
All_species_benchmark <- read.csv("All_Species_NSR_Pools_Timing_Benchmark.csv")
All_species_benchmark$transcriptome_size_MB=All_species_benchmark$transcriptome_size/1000000
PlotC <- ggplot(All_species_benchmark, aes(x=transcriptome_size_MB,y=time))+
  geom_point(color="red")+
  geom_smooth(method="lm",se=FALSE,color="black")+
  xlab("Transcriptome Size (MB)")+
  ylab("Run Time (s)")+
  theme_classic()
eqdata <- data.frame(x=All_species_benchmark$transcriptome_size_MB,y=All_species_benchmark$time)
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

PlotC_eq <- PlotC + geom_text(x = 2, y = 2000, label = lm_eqn(eqdata), parse = TRUE)
PlotC_eq
##########################3

SMU_hit_only <- read.csv("S_mutans_Cumulative_NSR_Pool_Genes_Hit_Only_5-28-21.csv")
SMU_hit_only$failed <- grepl("-",SMU_hit_only$randomer)
SMU_hit_only <- SMU_hit_only[which(SMU_hit_only$failed==FALSE), ]
new_run=1:nrow(SMU_hit_only)
SMU_hit_only$new_run <- new_run
PlotD <- ggplot(SMU_hit_only,aes(x=new_run,y=cumulative_genes_hit))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=1895,colour="black")+
  xlab("Pool Size")+
  ylab("Number of Genes Hit")+
  theme_classic()
PlotD

######################################################
SMU_log_only <- read.csv("S_mutans_Cumulative_NSR_Pool_Log_Term_6-1-21.csv")
SMU_log_only$failed <- grepl("-",SMU_log_only$randomer)
SMU_log_only <- SMU_log_only[which(SMU_log_only$failed==FALSE), ]
new_run=1:nrow(SMU_log_only)
SMU_log_only$new_run <- new_run
PlotE <- ggplot(SMU_log_only,aes(x=new_run,y=cumulative_hits))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=412260,colour="black")+
  xlab("Pool Size")+
  ylab("Total Hits")+
  theme_classic()
PlotE
cowplot::plot_grid(PlotD, PlotE, PlotA, align="h",nrow=1)

SMU_optimal <- read.csv("S_mutans_Cumulative_NSR_Pool_Optimal_Tuning_6-3-21.csv")
SMU_optimal$failed <- grepl("-",SMU_optimal$randomer)
SMU_optimal <- SMU_optimal[which(SMU_optimal$failed==FALSE), ]
new_run=1:nrow(SMU_optimal)
SMU_optimal$new_run <- new_run
optPlotA <- ggplot(SMU_optimal,aes(x=new_run,y=cumulative_genes_hit))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=1895,colour="black")+
  xlab("Pool Size")+
  ylab("Genes Hit")+
  theme_classic()
optPlotB <- ggplot(SMU_optimal,aes(x=new_run,y=cumulative_hits))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=412260,colour="black")+
  xlab("Pool Size")+
  ylab("Total Hits")+
  theme_classic()
optPlotC <- ggplot(SMU_optimal,aes(x=new_run,y=cumulative_dist_score))+
  geom_line(colour="red")+
  geom_vline(xintercept=453,colour="black")+
  geom_hline(yintercept=.337,colour="black")+
  xlab("Pool Size")+
  ylab("DUlq Score")+
  theme_classic()
cowplot::plot_grid(optPlotA,optPlotB, optPlotC, align="h",nrow=1)


################ E. coli Cumulative Plots for comparison with Brute Force #################
Ecoli_gene_only  <- read.csv("E_coli_Cumulative_NSR_Pool_Genes_Hit_Only_6-3-21.csv")
Ecoli_gene_only$failed <- grepl("-",Ecoli_gene_only$randomer)
Ecoli_gene_only <- Ecoli_gene_only[which(Ecoli_gene_only$failed==FALSE), ]
new_run=1:nrow(Ecoli_gene_only)
Ecoli_gene_only$new_run <- new_run
Ecoli_PlotA <- ggplot(Ecoli_gene_only,aes(x=new_run,y=cumulative_genes_hit))+
  geom_line(color="red")+
  geom_vline(xintercept = 437, color="black")+
  geom_hline(yintercept = 4436,color="black")+
  xlab("Pool Size")+
  ylab("Number of Genes Hit")+
  theme_classic()

Ecoli_uniformity_only <- read.csv("E_coli_Cumulative_NSR_Pool_Uniformity_Only_6-3-21.csv")
Ecoli_uniformity_only$failed <- grepl("-",Ecoli_uniformity_only$randomer)
Ecoli_uniformity_only <- Ecoli_uniformity_only[which(Ecoli_uniformity_only$failed==FALSE), ]
new_run=1:nrow(Ecoli_uniformity_only)
Ecoli_uniformity_only$new_run <- new_run
Ecoli_PlotB <- ggplot(Ecoli_uniformity_only,aes(x=new_run,y=cumulative_dist_score))+
  geom_line(color="red")+
  geom_vline(xintercept = 437, color="black")+
  geom_hline(yintercept = .299,color="black")+
  xlab("Pool Size")+
  ylab("DUlq Score")+
  theme_classic()
Ecoli_log_only <- read.csv("E_coli_Cumulative_NSR_Pool_Log_Term_Only_6-3-21.csv")
Ecoli_log_only$failed <- grepl("-",Ecoli_log_only$randomer)
Ecoli_log_only <- Ecoli_log_only[which(Ecoli_log_only$failed==FALSE), ]
new_run=1:nrow(Ecoli_log_only)
Ecoli_log_only$new_run <- new_run
Ecoli_PlotC <- ggplot(Ecoli_log_only,aes(x=new_run,y=cumulative_hits))+
  geom_line(color="red")+
  geom_vline(xintercept = 437, color="black")+
  geom_hline(yintercept = 719991,color="black")+
  xlab("Pool Size")+
  ylab("Total Hits")+
  theme_classic()

cowplot::plot_grid(Ecoli_PlotA, Ecoli_PlotC, Ecoli_PlotB, align="h",nrow=1)
plot_genome_timing_NSR <- function(file){
  data <- read.csv(file)
  data$transcriptome_sizeMB <- data$transcriptome_size/1e6
  data$timeMin <- data$time/60;
  PlotA <- ggplot2::ggplot(data,(aes(x=transcriptome_sizeMB,y=timeMin)))+
    geom_point(color="red")+
    xlab("Transcriptome Size (Mb)")+
    ylab("Time (min)")+
    geom_smooth(method="lm", se=FALSE,color="red")+
    theme_classic(base_size=12,base_family = "Arial")+
    ggrepel::geom_label_repel(aes(label = species),
                     box.padding   =1, 
                     point.padding = 0.1,
                     max.overlaps = 15,
                     segment.color = 'grey50')
  
  PlotA
}
plot_genome_timing_NSR("All_Species_NSR_Pools_Timing_Benchmark.csv")
plot_pool_size_timing_NSR <- function(file){
  data <- read.csv(file)
  data$cumulative_timeMin <- data$cumulative_time/60
  plotA <- ggplot(data,aes(x=pool_size,y=cumulative_timeMin))+
    geom_point(color="red")+
    theme_classic(base_size=12,base_family = "Arial")+
    xlab("Pool Size")+
    ylab("Time (min)")
  plotA
} 
plot_pool_size_timing_NSR("NSR_pool_size_timing_6_18_21.csv")
### Function to make the cumulative plots 
plot_cumulative_NSR_Experiment <- function(infile,outfile,BF_poolsize=453,BF_genes_hit=1895, BF_total_hits=412260, BF_interuniformity_score=.337,BF_intrauniformity_score=0.23){
  data <- read.csv(infile)
  data$failed <- grepl("-",data$randomer)
  data <- data[which(data$failed==FALSE), ]
  new_run=1:nrow(data)
  data$new_run <- new_run
  data$hits_thousands=data$cumulative_hits/1000;
  BF_total_hits_thousands=BF_total_hits/1000;
  im_width=5
  im_height=1.25
  font_size=6
  font="Arial"
  line_size=.25
  PlotA <- ggplot(data,aes(x=new_run,y=cumulative_genes_hit))+
    geom_line(color="red")+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_genes_hit,colour="black",linetype="dashed")+
    xlab("Pool Size")+
    ylab("Genes Hit")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotB <- ggplot(data,aes(x=new_run,y=hits_thousands))+
    geom_line(color="red")+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_total_hits_thousands,colour="black",linetype="dashed")+
    scale_y_continuous(labels = scales::comma)+
    xlab("Pool Size")+
    ylab("Total Hits (Thousands)")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotC <- ggplot(data,aes(x=new_run,y=cumulative_interuniformity_score))+
    geom_line(color="red")+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_interuniformity_score,colour="black",linetype="dashed")+
    xlab("Pool Size")+
    ylab("Intergene Uniformity")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotD <- ggplot(data,aes(x=new_run,y=cumulative_intrauniformity_score))+
    geom_line(color="red")+
    geom_vline(xintercept=BF_poolsize,color="black",linetype="dashed")+
    geom_hline(yintercept=BF_intrauniformity_score,linetype="dashed",colour="black")+
    xlab("Pool Size")+
    ylab("Intragene Uniformity")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
   cowplot::plot_grid(PlotA,PlotB,PlotC,PlotD, align="vh",nrow=1)


  ggsave(outfile,width=im_width,height=im_height,dpi=300,units="in")
}

plot_NSR_Performance_Data <- function(infile1,infile2,infile3,infile4,BF_poolsize=453,BF_genes_hit=1895, BF_total_hits=412260, BF_interuniformity_score=.337,BF_intrauniformity_score=0.23){
  ##Data1
  data1 <- read.csv(infile1)
  data1$failed <- grepl("-",data1$randomer)
  data1 <- data1[which(data1$failed==FALSE), ]
  new_run=1:nrow(data1)
  data1$new_run <- new_run
  data1$hits_thousands=data1$cumulative_hits/1000;
  ##Data2
  data2 <- read.csv(infile2)
  data2$failed <- grepl("-",data2$randomer)
  data2 <- data2[which(data2$failed==FALSE), ]
  new_run=1:nrow(data2)
  data2$new_run <- new_run
  data2$hits_thousands=data2$cumulative_hits/1000;
  ##Data1
  data3 <- read.csv(infile3)
  data3$failed <- grepl("-",data3$randomer)
  data3 <- data3[which(data3$failed==FALSE), ]
  new_run=1:nrow(data3)
  data3$new_run <- new_run
  data3$hits_thousands=data3$cumulative_hits/1000;
  ##Data1
  data4 <- read.csv(infile4)
  data4$failed <- grepl("-",data4$randomer)
  data4 <- data4[which(data4$failed==FALSE), ]
  new_run=1:nrow(data4)
  data4$new_run <- new_run
  data4$hits_thousands=data4$cumulative_hits/1000;
  BF_total_hits_thousands=BF_total_hits/1000;
  im_width=5
  im_height=1.25
  font_size=8
  font="Arial"
  line_size=.4
  red <- "#d73027"
  PlotA <- ggplot(data1,aes(x=new_run,y=cumulative_genes_hit))+
    geom_line(color=red)+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_genes_hit,colour="black",linetype="dashed")+
    xlab("Pool Size")+
    ylab("Genes Hit")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotB <- ggplot(data2,aes(x=new_run,y=hits_thousands))+
    geom_line(color=red)+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_total_hits_thousands,colour="black",linetype="dashed")+
    scale_y_continuous(labels = scales::comma)+
    xlab("Pool Size")+
    ylab("Total Hits (thousands)")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotC <- ggplot(data3,aes(x=new_run,y=cumulative_interuniformity_score))+
    geom_line(color=red)+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_interuniformity_score,colour="black",linetype="dashed")+
    xlab("Pool Size")+
    ylab("Intergene Uniformity")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotD <- ggplot(data4,aes(x=new_run,y=cumulative_intrauniformity_score))+
    geom_line(color=red)+
    geom_vline(xintercept=BF_poolsize,color="black",linetype="dashed")+
    geom_hline(yintercept=BF_intrauniformity_score,linetype="dashed",colour="black")+
    xlab("Pool Size")+
    ylab("Intragene Uniformity")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)

  return(cowplot::plot_grid(PlotA,PlotB,PlotC,PlotD, align="vh",nrow=1))
  
  

}

#genome_bench_file <- "All_Species_NSR_Pools_Timing_Benchmark.csv" 
#poolsize_bench_file <- "NSR_pool_size_timing_6_23_21.csv"
#plot_NSR_Benchmarking_Data(genome_bench_file,poolsize_bench_file)
plot_NSR_Benchmarking_Data <- function(infile1,infile2){
  font_size=8
  font="Arial"
  line_size=.4
  point_size=0.6
  red="#d73027"
  blue="#4575b4"
  data1 <- read.csv(infile1,stringsAsFactors = FALSE)
  data1$transcriptome_sizeMB <- data1$transcriptome_size/1e6
  data1$timeMin <- data1$time/60;
  data1$formatted_species <- paste0("italic('",str_replace_all(data1$species,"_",". "),"')")
  data1$label_species <- NA
  spec_label_idx <- c(1,2,3,4,21,10,18)
  data1$label_species[spec_label_idx] <- data1$formatted_species[spec_label_idx]
  PlotA <- ggplot2::ggplot(data1,(aes(x=transcriptome_sizeMB,y=timeMin,label=label_species)))+
    geom_point(color=red,size=point_size)+
    xlab("Transcriptome Size (Mb)")+
    ylab("Time (min)")+
    geom_smooth(method="lm", se=FALSE,color=red,size=point_size)+
    ggrepel::geom_text_repel(
      segment.size=0.4,
      size=8/3,
      box.padding   =.5, 
      point.padding = 0.2,
      max.overlaps = 10,
      segment.color = 'grey50',
      parse=TRUE,
      na.rm=TRUE,
      force=15,
      min.segment.length = 0,
      nudge_x = c(-1,-1,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0),
      nudge_y = c(2,3,-1,-2,0,0,0,0,0,2,0,0,0,0,0,0,0,-2,0,0,-2,0,0,0,0))+
    theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)

  data2 <- read.csv(infile2)
  data2$without_intraMin <- data2$without_intra/60
  data2$with_intraMin <- data2$with_intra/60
  data2 <- data2[seq(1,nrow(data2),by=20), ]
  PlotB <- ggplot(data2)+
    geom_point(aes(x=run,y=without_intraMin),color=blue,size=point_size)+
    geom_smooth(aes(x=run,y=without_intraMin),method="lm",se=FALSE,color=blue,size=point_size)+
    geom_point(aes(x=run,y=with_intraMin),color="black",size=point_size)+
    geom_smooth(aes(x=run,y=with_intraMin),method="lm",se=FALSE,color="black",size=point_size)+
    theme_classic(base_size=font_size,base_family = font,base_line_size = 1.1*line_size)+
    xlab("Pool Size")+
    ylab("Time (min)")
  plts <- cowplot::plot_grid(PlotA,PlotB,NULL, rel_widths = c(2,1,1),align="vh",nrow=1)
  return(plts)
}

plot_NSR_figure <- function(outfile){
  scale_factor=6
  im_width=scale_factor
  im_height=scale_factor/2
  
  genome_bench_file <- "All_Species_NSR_Pools_Timing_Benchmark.csv" 
  poolsize_bench_file <- "NSR_pool_size_timing_6_23_21.csv"
  infile1 <-"S_mutans_Cumulative_NSR_CCD_Run_12_6_30_21.csv"
  infile2 <-"S_mutans_Cumulative_NSR_CCD_Run_12_6_30_21.csv"
  infile3 <-"S_mutans_Cumulative_NSR_CCD_Run_12_6_30_21.csv"
  infile4 <-"S_mutans_Cumulative_NSR_CCD_Run_12_6_30_21.csv"
  row1 <- plot_NSR_Performance_Data(infile1,infile2,infile3,infile4)
  row2 <- plot_NSR_Benchmarking_Data(genome_bench_file,poolsize_bench_file)
  cowplot::plot_grid(row1,row2,nrow=2,ncol=1,rel_widths=c(1,1),align="h")
  ggsave(outfile,width=im_width,height=im_height,dpi=300,units="in")
}

plot_NSR_figure("test.tiff")

#plot_cumulative_NSR("S_mutans_Cumulative_NSR_Pool_Best_Sim_Option_6_8_21.csv")
#plot_cumulative_NSR("S_mutans_Cumulative_NSR_Pool_Best_Sim_Option_6_8_21_Gene.csv")
#plot_cumulative_NSR("S_mutans_Cumulative_NSR_Pool_Best_Sim_Option_6_8_21_Uniformity.csv")

plot_cumulative_NSR_Experiment("S_mutans_Cumulative_NSR_Pool_All_Rewards_6_22_21_1.csv","Cumulative_NSR_Pools_6_22_21_1.tiff")
plot_cumulative_NSR_Experiment("S_mutans_Cumulative_NSR_Pool_All_Rewards_6_22_21_2.csv","Cumulative_NSR_Pools_6_22_21_2.tiff")
plot_cumulative_NSR_Experiment("S_mutans_Cumulative_NSR_Pool_All_Rewards_6_22_21_3.csv","Cumulative_NSR_Pools_6_22_21_3.tiff")
plot_cumulative_NSR_Experiment("S_mutans_Cumulative_NSR_Pool_All_Rewards_6_22_21_4.csv","Cumulative_NSR_Pools_6_22_21_4.tiff")
