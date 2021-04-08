library(ggplot2)
fsts_sim<- do.call("rbind",list(low_fsts,marg_fsts,high_fsts,extreme_fsts))
fsts_sim$Scenario <- word(fsts_sim$FileName,1,sep="_")
fsts_sim$Scen_F = factor(fsts_sim$Scenario, levels=c('low','marg','high','extreme'))
fsts_sim$Lab<-c(rep("A",440),rep("B",440),rep("C",440),rep("D",440))
struc_sim<- do.call("rbind",list(low_struc_mats,marg_struc_mats,high_struc_mats,extreme_struc_mats))
struc_sim$Scenario <- word(struc_sim$files,2,sep="_")
struc_sim$Scen_F = factor(struc_sim$Scenario, levels=c('low','marg','high','extreme'))

pcst_sim<- do.call("rbind",list(pcst_lowPS,pcst_margPS,pcst_highPS,pcst_extremePS))
pcst_sim$Scenario <- word(pcst_sim$FileNames,1,sep="_")
pcst_sim$Scen_F = factor(pcst_sim$Scenario, levels=c('low','marg','high','extreme'))
pcst_sim$Lab<-c(rep("A",440),rep("B",440),rep("C",440),rep("D",440))
library(rstatix);library(ggpubr)
(pcst_sim %>% filter(Scen_F == "marg") %>% pairwise_wilcox_test(
  Fst ~ Filter, paired = TRUE,
  p.adjust.method = "bonferroni"
))

View(struc_sim %>% filter(Scen_F == "extreme") %>% pairwise_wilcox_test(
  struc_mats ~ filter, paired = TRUE,
  p.adjust.method = "bonferroni"
))

ggplot(fsts_sim,aes(x=InferredFst,y=Filt))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter")

ggplot(struc_sim,aes(x=struc_mats,y=filter))+  #theme_minimal() 
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Nucleotide Distance") 
ggplot(pcst_sim,aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("PCst")

seal_fsts$DatType<-"seal"
isopod_fsts$DatType<-"isopod"
zebra_fsts$DatType<-"zebra"

real_fst<-rbind(seal_fsts,isopod_fsts,zebra_fsts)

real_fst$Scen_F = factor(real_fst$DatType, levels=c('seal','zebra','isopod'))

ggplot(real_fst,aes(x=Fst,y=Filter))+  #
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter")


real_pcst<-rbind(pcst_sealPS,pcst_isopodPS,pcst_zebraPS)
real_pcst$DatType<-c(rep("seal",40),rep("isopod",40),rep("zebra",40))
real_pcst$Scen_F = factor(real_pcst$DatType, levels=c('seal','zebra','isopod'))

ggplot(real_pcst,aes(x=V1,y=Filter))+  #
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red", bandwidth=0.006) +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") +xlab ("PCst")
# Real structure
struc_mats_seal$DatType<-"seal"
isopod_struc_mats$DatType<-"isopod"
struc_mats_zebra$DatType<-"zebra"
struc_dat_real<-rbind(struc_mats_seal,isopod_struc_mats,struc_mats_zebra)
struc_dat_real$Scen_F = factor(struc_dat_real$DatType, levels=c('seal','zebra','isopod'))
ggplot(struc_dat_real,aes(x=struc_mats,y=filter))+  #theme_minimal() 
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.0000000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Nucleotide Distance")


ggplot(isopod_polyploidy_struc_mats,aes(x=struc_mats,y=filter))+  #theme_minimal() 
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.0000000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red",bandwidth=0.005) +
  theme_bw() + ylab("HWE Filter") + xlab("Nucleotide Distance") 



#### RANDOM
random_dat<-rbind(pcst_low_randPS,pcst_extreme_randPS)
random_dat$Scen_F <- c(rep("low",440),rep("extreme",440))
random_dat$Scen_F = factor(random_dat$Scen_F, levels=c('low','extreme'))
ggplot(random_dat,aes(x=V1,y=Filter))+  #
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("PCst")

ggplot(random_dat,aes(x=Fst,y=Filter))+  #
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Fst")

