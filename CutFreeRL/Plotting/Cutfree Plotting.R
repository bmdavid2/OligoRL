setwd("~/Documents/GitHub/OligoRL/CutFreeRL/Plotting")
library(ggplot2)
library(dplyr)

results_sites <- read.csv("./rollout_results_sites_4_1_21.csv")
results_length <- read.csv("./rollout_results_length_4_15_21.csv")
interaction_data <- read.csv("./Cutfree_interaction_data.csv")
#results_sites <- results_sites[order(results_sites$n_sites), ]
#resutls_sites <- results_sites[results_sites$n_sites<=10]
#results_sites <- results_sites[1:80,]
plot_cutfree_data <- function() {
  font_size=8
  font="Arial"
  line_size=.4
  point_size=0.6
  red="#d73027"
  blue="#4575b4"
  results_sites <- read.csv("./rollout_results_sites_4_1_21.csv")
  results_length <- read.csv("./rollout_results_length_4_15_21.csv")
  interaction_data <- read.csv("./Cutfree_interaction_data.csv")
  plot1 <- ggplot2::ggplot(data=results_sites) +
    ggplot2::scale_y_log10(breaks=c(0.01,0.1,1,10,100),labels=c("0.01","0.1","1","10","100")) +
    ggplot2::geom_point(aes(x=n_sites,y=time),color="black",size=point_size) +
    ggplot2::geom_point(aes(x=n_sites,y=random_runtime),color=red,size=point_size)+
    ggplot2::xlab("Restriction Sites") +
    ggplot2::ylab("Run Time (s)") +
    ggplot2::geom_smooth(aes(x=n_sites,y=time),method="lm", se=FALSE,color="black")+
    ggplot2::geom_smooth(aes(x=n_sites,y=random_runtime),method="lm", se=FALSE,color=red)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)
  
  plot2 <- ggplot2::ggplot(data=results_sites) +
    ggplot2::scale_y_log10(breaks=c(1e9,1e10,1e11),labels=c(expression(paste("10"^"9")),expression(paste("10"^"10")),expression(paste("10"^"11")))) +
    ggplot2::geom_point(aes(x=n_sites,y=degeneracy),color="black",size=point_size) +
    ggplot2::geom_point(aes(x=n_sites,y=random_objval),color=red,size=point_size)+
    ggplot2::xlab("Restriction Sites") +
    ggplot2::ylab("Possible Sequences") +
    ggplot2::geom_smooth(aes(x=n_sites,y=degeneracy),method="lm", se=FALSE,color="black")+
    ggplot2::geom_smooth(aes(x=n_sites,y=random_objval),method="lm", se=FALSE,color=red)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)
  
  plot0 <- ggplot2::ggplot(data=results_length) +
    ggplot2::scale_y_log10(breaks=c(0.1,1,10,100),labels=c("0.1","1","10","100")) +
    ggplot2::geom_point(aes(x=oligo_lengths, y=time),color="black",size=point_size) +
    ggplot2::geom_point(aes(x=oligo_lengths, y=random_runtime),color=red,size=point_size) +
    ggplot2::xlab("Randomer Length (bp)") +
    ggplot2::ylab("Run Time (s)") +
    ggplot2::geom_smooth(aes(x=oligo_lengths,y=time),method="lm", se=FALSE,color="black")+
    ggplot2::geom_smooth(aes(x=oligo_lengths,y=random_runtime),method="lm", se=FALSE,color=red)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)
  
  #results %>%
    #dplyr::group_by(n_sites) %>%
    #dplyr::select(n_sites, time, degeneracy) %>%
    #dplyr::summarise_all(dplyr::funs(mean, sd))
  
  cowplot::plot_grid(plot1, plot0,plot2,NULL,rel_widths = c(1,1,1,0.75), align= "h", nrow=1 )
  scaling_factor=5.5
  ggsave("CutfreeRL_Figure.tiff",width=scaling_factor,height=scaling_factor/3.75,dpi=300,units="in")
}
plot_cutfree_data()
## Test cutfree runtime significance 
summary(logcutfree.model <- lm(log(time) ~ n_sites,data=results_sites))
MASS::boxcox(time ~ n_sites, data=results_sites)
plot(results_sites$n_sites,log(results_sites$time))

## Test RL runtime significance 
summary(rlcutfree.model <- lm((random_runtime) ~ n_sites,data=results_sites))
MASS::boxcox(random_runtime ~ n_sites, data=results_sites)
plot(results_sites$n_sites,results_sites$random_runtime)

## test cutfree solution decay 
summary(solution.cutfree <- lm(log(degeneracy)~n_sites,data=results_sites))
MASS::boxcox(degeneracy ~ n_sites, data=results_sites)
## Test RL solution decay 
summary(solution.RL <- lm(log(random_objval)~n_sites,data=results_sites))
MASS::boxcox(random_objval ~ n_sites, data=results_sites)
##Test interaction between nsites and algorithm
summary(solution.interaction <- lm(log(degeneracy)~n_sites+algorithm+n_sites*algorithm,data=interaction_data))
plot(model)
plot(log(interaction_data$degeneracy),predict(solution.interaction))

##test cutfree oligo length 
summary(length.cutfree <- lm(log(time) ~ oligo_lengths,data=results_length))
summary(length.RL <- lm(log(random_runtime) ~oligo_lengths,data=results_length))
MASS::boxcox(random_runtime ~ oligo_lengths, data=results_length)
MASS::boxcox(time ~oligo_lengths,data=results_length)
