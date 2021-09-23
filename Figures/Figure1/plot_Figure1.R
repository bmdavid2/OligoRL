
library(ggplot2)
library(dplyr)

results_sites <- read.csv("./rollout_results_sites.csv")
results_length <- read.csv("./rollout_results_length.csv")
interaction_data <- read.csv("./Cutfree_interaction_data.csv")

plot_cutfree_data <- function(results_sites,results_length) {
  font_size=8
  font="Arial"
  line_size=.4
  point_size=0.6
  red="#d73027"
  blue="#4575b4"
  plot1 <- ggplot2::ggplot(data=results_sites) +
    ggplot2::scale_y_log10(breaks=c(0.01,0.1,1,10,100),labels=c("0.01","0.1","1","10","100")) +
    ggplot2::geom_point(aes(x=n_sites,y=time),color="black",size=point_size) +
    ggplot2::geom_point(aes(x=n_sites,y=random_runtime),color=red,size=point_size)+
    ggplot2::xlab("Restriction Sites") +
    ggplot2::ylab("Runtime (s)") +
    ggplot2::labs(color="black")+
    ggplot2::geom_smooth(aes(x=n_sites,y=time),method="lm", se=FALSE,color="black")+
    ggplot2::geom_smooth(aes(x=n_sites,y=random_runtime),method="lm", se=FALSE,color=red)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)
  
  plot2 <- ggplot2::ggplot(data=results_sites) +
    ggplot2::scale_y_log10(breaks=c(1e9,1e10,1e11),labels=c(expression(paste("10"^"9")),expression(paste("10"^"10")),expression(paste("10"^"11")))) +
    ggplot2::geom_point(aes(x=n_sites,y=degeneracy),color="black",size=point_size) +
    ggplot2::geom_point(aes(x=n_sites,y=random_objval),color=red,size=point_size)+
    ggplot2::xlab("Restriction Sites") +
    ggplot2::ylab("Possible Sequences") +
    ggplot2::labs(color="black")+
    ggplot2::geom_smooth(aes(x=n_sites,y=degeneracy),method="lm", se=FALSE,color="black")+
    ggplot2::geom_smooth(aes(x=n_sites,y=random_objval),method="lm", se=FALSE,color=red)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)
  
  plot0 <- ggplot2::ggplot(data=results_length) +
    ggplot2::scale_y_log10(breaks=c(0.1,1,10,100),labels=c("0.1","1","10","100")) +
    ggplot2::geom_point(aes(x=oligo_lengths, y=time),color="black",size=point_size) +
    ggplot2::geom_point(aes(x=oligo_lengths, y=random_runtime),color=red,size=point_size) +
    ggplot2::xlab("Randomer Length (bp)") +
    ggplot2::ylab("Runtime (s)") +
    ggplot2::geom_smooth(aes(x=oligo_lengths,y=time),method="lm", se=FALSE,color="black")+
    ggplot2::geom_smooth(aes(x=oligo_lengths,y=random_runtime),method="lm", se=FALSE,color=red)+
    ggplot2::theme_classic(base_size=font_size,base_family = font,base_line_size = line_size)

  
  #results %>%
    #dplyr::group_by(n_sites) %>%
    #dplyr::select(n_sites, time, degeneracy) %>%
    #dplyr::summarise_all(dplyr::funs(mean, sd))
  
  cowplot::plot_grid(plot1, plot0,plot2,NULL,rel_widths = c(1,1,1,0.75), align= "h", nrow=1 )
  scaling_factor=5.5
  ggsave("./Figure1.tiff",width=scaling_factor,height=scaling_factor/3.75,dpi=300,units="in")
}
plot_cutfree_data(results_sites,results_length)

##Test interaction between nsites and algorithm
summary(solution.interaction <- lm(log(degeneracy)~n_sites+algorithm+n_sites*algorithm,data=interaction_data))
plot(solution.interaction)
plot(log(interaction_data$degeneracy),predict(solution.interaction))


