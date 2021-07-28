setwd("~/Documents/GitHub/OligoRL/NSR_RL/Plotting")
library(ggplot2)
library(stringr)


plot_NSR_Performance_Data <- function(infile1,infile2,infile3,infile4,BF_poolsize=453,BF_pct_genes_hit=100, BF_total_hits=412260, BF_interuniformity_score=.337,BF_intrauniformity_score=0.237){
  ##Data1
  data1 <- read.csv(infile1)
  data1$failed <- grepl("-",data1$randomer)
  data1 <- data1[which(data1$failed==FALSE), ]
  new_run=1:nrow(data1)
  data1$new_run <- new_run
  data1$hits_thousands=data1$cumulative_hits/1000;
  data1$hits_per_gene=data1$cumulative_hits/1895;
  data1$pct_genes_hit=data1$cumulative_genes_hit/1895*100;
  ##Data2
  data2 <- read.csv(infile2)
  data2$failed <- grepl("-",data2$randomer)
  data2 <- data2[which(data2$failed==FALSE), ]
  new_run=1:nrow(data2)
  data2$new_run <- new_run
  data2$hits_thousands=data2$cumulative_hits/1000;
  data2$hits_per_gene=data2$cumulative_hits/1895;
  data2$pct_genes_hit=data2$cumulative_genes_hit/1895*100;
  ##Data1
  data3 <- read.csv(infile3)
  data3$failed <- grepl("-",data3$randomer)
  data3 <- data3[which(data3$failed==FALSE), ]
  new_run=1:nrow(data3)
  data3$new_run <- new_run
  data3$hits_thousands=data3$cumulative_hits/1000;
  data3$hits_per_gene=data3$cumulative_hits/1895;
  data3$pct_genes_hit=data3$cumulative_genes_hit/1895*100;
  ##Data1
  data4 <- read.csv(infile4)
  data4$failed <- grepl("-",data4$randomer)
  data4 <- data4[which(data4$failed==FALSE), ]
  new_run=1:nrow(data4)
  data4$new_run <- new_run
  data4$hits_thousands=data4$cumulative_hits/1000;
  data4$hits_per_gene=data4$cumulative_hits/1895;
  data4$pct_genes_hit=data4$cumulative_genes_hit/1895*100;
  BF_total_hits_thousands=BF_total_hits/1000;
  BF_avg_hits_per_gene=BF_total_hits/1895
  im_width=5
  im_height=1.25
  font_size=8
  font="Arial"
  line_size=.4
  red <- "#d73027"
  PlotA <- ggplot(data1,aes(x=new_run,y=pct_genes_hit))+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_pct_genes_hit,colour="black",linetype="dashed")+
    geom_line(color=red)+
    xlab("Pool Size")+
    ylab("Genes Hit (%)")+
    scale_y_continuous(labels= scales::number_format(accuracy=1),limits=c(86,100))+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotB <- ggplot(data2,aes(x=new_run,y=hits_per_gene))+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_avg_hits_per_gene,colour="black",linetype="dashed")+
    geom_line(color=red)+
    scale_y_continuous(labels = scales::comma)+
    xlab("Pool Size")+
    ylab("Average Hits/Gene")+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotC <- ggplot(data3,aes(x=new_run,y=cumulative_interuniformity_score))+
    geom_vline(xintercept=BF_poolsize,colour="black",linetype="dashed")+
    geom_hline(yintercept=BF_interuniformity_score,colour="black",linetype="dashed")+
    geom_line(color=red)+
    xlab("Pool Size")+
    ylab("Intergene Uniformity")+
    scale_y_continuous(labels=c("0.28","0.30","0.32","0.34","0.36"),limits=c(0.26,0.36),breaks=c(0.28,0.30,0.32,0.34,0.36))+
    theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  PlotD <- ggplot(data4,aes(x=new_run,y=cumulative_intrauniformity_score))+
    geom_vline(xintercept=BF_poolsize,color="black",linetype="dashed")+
    geom_hline(yintercept=BF_intrauniformity_score,linetype="dashed",colour="black")+
    geom_line(color=red)+
    xlab("Pool Size")+
    ylab("Intragene Uniformity")+
    scale_y_continuous(labels=c("0.16","0.18","0.20","0.22","0.24"),limits=c(0.14,0.24),breaks=c(0.16,0.18,0.20,0.22,0.24))+
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
    geom_smooth(method="lm", se=FALSE,color=red,size=point_size)+
    geom_point(color=red,size=point_size)+
    xlab("Transcriptome Size (Mb)")+
    ylab("Runtime (min)")+
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
    geom_smooth(aes(x=run,y=without_intraMin),method="lm",se=FALSE,color=blue,size=point_size)+
    geom_point(aes(x=run,y=without_intraMin),color=blue,size=point_size)+
    geom_smooth(aes(x=run,y=with_intraMin),method="lm",se=FALSE,color="black",size=point_size)+
    geom_point(aes(x=run,y=with_intraMin),color="black",size=point_size)+
    theme_classic(base_size=font_size,base_family = font,base_line_size = 1.1*line_size)+
    xlab("Pool Size")+
    ylab("Runtime (min)")
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

plot_NSR_figure("NSR_figure.tiff")

