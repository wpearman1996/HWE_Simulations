setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/k12_extreme/")
k12_files<-list.files(pattern="*out_f")
#k12_files<-k12_files[grepl("k12",k12_files)]
k12_files<-k12_files[!grepl("4817126967349",k12_files)]


library(tidyverse)
variance_ln<-function(file){
  z<-readLines(file)
  z<-z[grepl("Mean value of ln likelihood",z)]
  as.numeric(word(z,2,sep="= "))
}
variances<-lapply(k12_files,variance_ln)
names(variances)<-k12_files
variances<-do.call("rbind",variances)
variances<-data.frame(variances,k12_files)
variances$filter<-ifelse(grepl("nohwe",variances$k12_files),"No Filter",
                         ifelse(grepl("out_any",variances$k12_files),"Out Any",
                                ifelse(grepl("out_across",variances$k12_files),
                                       "Out Combo", ifelse(grepl("within",variances$k12_files), 
                                                           "Out Within", "Out All"))))
variances$Scenario<-ifelse(grepl("high",variances$k12_files),"High",
                           ifelse(grepl("extreme",variances$k12_files),"Extreme",
                                  ifelse(grepl("marg",variances$k12_files), "Marginal", 
                                         ifelse(grepl("panmictic",variances$k12_files),"Panmictic","Low")))) 
variances$Linkage<-ifelse(grepl("lowlink",variances$k12_files),"Low Linkage","High Linkage")
variances$Scenario = factor(variances$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))                      
variances$Seed<-gsub('.*_([0-9]+).*','\\1',variances$k12_files)
variances<-variances[!grepl("*gz",variances$k12_files),]
variances$VarStandVariable<-paste(variances$Seed,variances$Linkage,variances$filter,variances$Scenario,sep="_")

stats::aggregate(.~variances$VarStandVariable,variances$variances,mean)
maxed<-aggregate(variances$variances, list(variances$VarStandVariable), max)
mined<-aggregate(variances$variances, list(variances$VarStandVariable), min)
maxed$Dif<-maxed$x-mined$x
library(ggridges)
variances$PerGroup_MaxVar<-maxed$Dif[match(variances$VarStandVariable,maxed$Group.1)]
variances$StandVar<-variances$variances/variances$PerGroup_MaxVar
struc_variance_plot<- variances %>% 
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=variances,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Linkage,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Mean log likelihood")
struc_variance_plot




k12_dat<-lapply(k12_files,read_struc_nucmat_k12)

k12_dat<-do.call("rbind",k12_dat)
k12_dat<-data.frame(k12_dat,k12_files)
k12_dat$filter<-ifelse(grepl("nohwe",k12_dat$k12_files),"No Filter",
                      ifelse(grepl("out_any",k12_dat$k12_files),"Out Any",
                             ifelse(grepl("out_across",k12_dat$k12_files),
                                    "Out Combo", ifelse(grepl("within",k12_dat$k12_files), 
                                                        "Out Within", "Out All"))))
k12_dat$Scenario<-ifelse(grepl("high",k12_dat$k12_files),"High",
                        ifelse(grepl("extreme",k12_dat$k12_files),"Extreme",
                               ifelse(grepl("marg",k12_dat$k12_files), "Marginal", 
                                      ifelse(grepl("panmictic",k12_dat$k12_files),"Panmictic","Low")))) 
k12_dat$Linkage<-ifelse(grepl("lowlink",k12_dat$k12_files),"Low Linkage","High Linkage")
k12_dat$Scenario = factor(k12_dat$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))                      
k12_dat$Seed<-gsub('.*_([0-9]+).*','\\1',k12_dat$k12_files)
k12_dat<-k12_dat[!grepl(".gz",k12_dat$k12_files),]



Extreme_k12_struc<- k12_dat %>% #filter(Linkage == "Low Linkage") %>%
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=k12_dat,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Linkage,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Average Nucleotide Distance")
Extreme_k12_struc

setwd("../Figures/")
ggsave(struc_variance_plot, filename = "Extreme_k12_strucvariance.pdf",width = 6,height=3.86)

ggsave(Extreme_k12_struc, filename = "Extreme_k12_struc.pdf",width = 6,height=3.86)
