plot_grid(low_fst_means + ylab("Filter") + scale_y_discrete(labels=c("No HWE" = "No Filter"))
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          marg_fst_means + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank()) 
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          high_fst_means + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank())
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          extreme_fst_means + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank())
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          nrow = 1,labels=c("A","B","C","D"))

plot_grid(low_pcst_all + ylab("Filter") + scale_y_discrete(labels=c("No HWE" = "No Filter"))
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          marg_pcst_all + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank()) 
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          high_pcst_all + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank())
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          extreme_pcst_all + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank())
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          nrow = 1,labels=c("A","B","C","D"))

plot_grid(low_pcst_noacross + ylab("Filter") + scale_y_discrete(labels=c("No HWE" = "No Filter"))
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          marg_pcst_noacross + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank()) 
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          high_pcst_noacross + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank())
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          extreme_pcst_noacross + theme(axis.text.y=element_blank()) + theme(axis.title.y = element_blank())
          + geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x)),
          nrow = 1,labels=c("A","B","C","D"))
library("ggpubr")
ggpubr::ggarrange(low_fst_means,marg_fst_means,high_fst_means,extreme_fst_means, # list of plots
                  labels = "AUTO", # labels
                  common.legend = T, # COMMON LEGEND
                  legend = "bottom", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  ncol = 4)  # number of rows

                  