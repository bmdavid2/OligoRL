setwd("~/Documents/GitHub/OligoRL/Oligo_Compressor/Plotting")
library(ggplot2)
library(readxl)
library(dplyr)


plotoligocompressor <-  function() {
  font_size=8
  font="Arial"
  line_size=.4
  point_size=0.6
  red="#d73027"
  blue="#4575b4"
  poolsizedata <- read.csv('oligo_compressor_benchmark_5_11_21.csv')
  poolsizedata$timeMin=poolsizedata$time/60
  recompdata <- read.csv('recompression_experiment_4-23-21_nsims100.csv')
  lengthsdata <- read.csv('oligo_compressor_benchmark_lengths_5_12_21.csv')
  recompdata <- recompdata[order(recompdata$decompressedsize), ]
  recompdata$order <- c(1:nrow(recompdata))
  recompdata <- recompdata[1:(nrow(recompdata)-7), ]
  plot1 <- ggplot(poolsizedata)+
    ggplot2::geom_point(aes(x=poolsize,y=compressedsize),color=red,size=point_size) +
    ggplot2::geom_point(aes(x=poolsize,y=poolsize),color=blue,size=point_size)+
    ggplot2::xlab("Original Pool Size") +
    ggplot2::ylab("Pool Size") +
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  plot2 <- ggplot(recompdata)+
    ggplot2::geom_point(aes(x=order,y=recompressedsize),color=red,size=point_size)+
    ggplot2::geom_point(aes(x=order,y=poolsize), color="black",size=point_size)+
    ggplot2::geom_point(aes(x=order,y=decompressedsize),color=blue,size=point_size)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)+
    ggplot2::xlab("Run Number") +
    ggplot2::ylab("Pool Size") 
  
  plot3 <-  ggplot(poolsizedata)+
    ggplot2::geom_point(aes(x=poolsize,y=timeMin),color=red,size=point_size)+
    ggplot2::xlab("Original Pool Size") +
    ggplot2::ylab("Run Time (s)") +
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)
  #plot4 <- ggplot(lengthsdata)+
    #ggplot2::geom_point(aes(x=oligo_length,y=time),color="red")+
    #ggplot2::xlab("Oligo Length")+
    #ggplot2::ylab("Runtime [s]")+
    #ggplot2::theme_classic()
  
  
  cowplot::plot_grid(plot1,plot3,plot2,NULL,rel_widths = c(1,1,1,.75),nrow=1,align="vh")
  ggsave("Oligo_Compressor_Figure_Plots.tiff",width=5.5,height=5.5/3.75,dpi=300,units="in")
}
plotoligocompressor()

plotrecompression <- function(){
  red="#d73027"
  blue="#4575b4"
  recomp1000 <- read.csv('recompression_3_19_21_add_nsims100.csv')
  recomp1000 <- recomp1000[order(recomp1000$decompressedsize), ]
  recomp1000$order <- c(1:nrow(recomp1000))
  recomp1000 <- recomp1000[1:(nrow(recomp1000)-4), ]
  recomp100 <- read.csv('recompression_experiment_4-23-21_nsims100.csv')
  recomp100 <- recomp100[order(recomp100$decompressedsize), ]
  recomp100$order <- c(1:nrow(recomp100))
  recomp100 <- recomp100[1:(nrow(recomp100)-7), ]
  plot1 <- ggplot(recomp1000)+
    ggplot2::geom_point(aes(x=order,y=recompressedsize),color=red)+
    ggplot2::geom_point(aes(x=order,y=poolsize),color="black")+
    ggplot2::geom_point(aes(x=order,y=decompressedsize),color=blue)+
    ggplot2::theme_classic()+
    ggplot2::xlab("Run Number (nsims=1000)")+
    ggplot2::ylab("Pool Size")
  plot2 <- ggplot(recomp1000)+
    ggplot2::geom_point(aes(x=order,y=recompressed100),color=red)+
    ggplot2::geom_point(aes(x=order,y=poolsize),color="black")+
    ggplot2::geom_point(aes(x=order,y=decompressedsize),color=blue)+
    ggplot2::theme_classic()+
    ggplot2::xlab("Run Number (nsims=100)")+
    ggplot2::ylab("Pool Size")
  cowplot::plot_grid(plot1,plot2, align="h",nrow=1)
  
}
plotrecompression()

