
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
  poolsizedata <- read.csv('oligo_compressor_benchmark.csv')
  poolsizedata$timeMin=poolsizedata$time/60
  recompdata <- read.csv('recompression_experiments.csv')
  recompdata <- recompdata[order(recompdata$decompressedsize), ]
  recompdata$order <- c(1:nrow(recompdata))
  recompdata$efficiency <- (recompdata$decompressedsize-recompdata$recompressedsize)/(recompdata$decompressedsize-recompdata$poolsize)
  meaneff <- mean(recompdata$efficiency)
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
  print("Mean Efficiency")
  meaneff
  cowplot::plot_grid(plot1,plot3,plot2,NULL,rel_widths = c(1,1,1,.75),nrow=1,align="vh")
  ggsave("Figure2.tiff",width=5.5,height=5.5/3.75,dpi=300,units="in")
}


plotrecompression <- function(){
  font_size=8
  font="Arial"
  line_size=.4
  point_size=0.6
  red="#d73027"
  blue="#4575b4"
  recomp1000 <- read.csv('recompression_experiments.csv')
  recomp1000 <- recomp1000[order(recomp1000$decompressedsize), ]
  recomp1000$order <- c(1:nrow(recomp1000))
  recomp1000 <- recomp1000[1:(nrow(recomp1000)), ]
  plot1 <- ggplot(recomp1000)+
    ggplot2::geom_point(aes(x=order,y=recompressedsize),color=red,size=point_size)+
    ggplot2::geom_point(aes(x=order,y=poolsize),color="black",size=point_size)+
    ggplot2::geom_point(aes(x=order,y=decompressedsize),color=blue,size=point_size)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)+
    ggplot2::xlab("Run Number")+
    ggplot2::ylab("Pool Size")+
    ggplot2::labs(subtitle="nsims=1000")
  plot2 <- ggplot(recomp1000)+
    ggplot2::geom_point(aes(x=order,y=recompressed100),color=red,size=point_size)+
    ggplot2::geom_point(aes(x=order,y=poolsize),color="black",size=point_size)+
    ggplot2::geom_point(aes(x=order,y=decompressedsize),color=blue,size=point_size)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size=line_size)+
    ggplot2::xlab("Run Number")+
    ggplot2::ylab("Pool Size")+
    ggplot2::labs(subtitle="nsims=100")
  cowplot::plot_grid(plot1,plot2, align="h",nrow=1)
  ggsave("nsims_comparison.tiff",width=3,height=1.6,dpi=300,units="in")
}
plotoligocompressor()
plotrecompression()

