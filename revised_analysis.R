setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/pca_structure_output/")
csvs <- list.files(path = "./",pattern="*.csv",
                           recursive = T,full.names = T)
csvs<-csvs[grepl("subrep",csvs)]
csvs_pca<-csvs[grepl("*PCA*",csvs)]
csvs_fst<-csvs[grepl("*Fst*",csvs)]

pcas<-lapply(csvs_pca,read.csv)
names(pcas)<-csvs_pca
pop_list<-c(rep("p1",30),rep("p2",30),rep("p3",30),rep("p4",30),rep("p5",30),rep("p6",30))
pcst_vals <- lapply(pcas,pc_st,pop_list)
pcst_vals<-do.call("rbind",pcst_vals)
pcst_vals<-as.data.frame(pcst_vals)
pcst_vals$File<-rownames(pcst_vals)
pcst_vals$Filter<-ifelse(grepl("nohwe",pcst_vals$File),"No Filter",
                         ifelse(grepl("out_any",pcst_vals$File),"Out Any",
                                ifelse(grepl("out_across",pcst_vals$File),
                                       "Out Combo", ifelse(grepl("within",pcst_vals$File),"Out Within","Out All")))) 
pcst_vals$Scenario<-ifelse(grepl("high",pcst_vals$File),"High",
                         ifelse(grepl("extreme",pcst_vals$File),"Extreme",
                                ifelse(grepl("marg",pcst_vals$File), "Marginal", 
                                       ifelse(grepl("panmictic",pcst_vals$File),"Panmictic","Low")))) 
pcst_vals$Linkage<-ifelse(grepl("lowlink",pcst_vals$File),"Low Linkage","High Linkage")
pcst_vals<-pcst_vals[!grepl("4817126967349",pcst_vals$File),]
pcst_vals$Scenario = factor(pcst_vals$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))
library(ggplot2)


library(rstatix)            
View(pcst_vals %>% filter(Scenario == "Extreme") %>%
       filter(Linkage == "Low Linkage") %>% pairwise_wilcox_test(
  V1 ~ Filter, paired = TRUE,
  p.adjust.method = "BH"
))
            

fsts<-lapply(csvs_fst,read.csv)
names(fsts)<-csvs_fst
fst_vals <- lapply(fsts,function(x){mean(as.matrix(x[2:7]),na.rm=T)})
fst_vals<-do.call("rbind",fst_vals)
fst_vals<-as.data.frame(fst_vals)
fst_vals$File<-rownames(fst_vals)
fst_vals$Filter<-ifelse(grepl("nohwe",fst_vals$File),"No Filter",
                         ifelse(grepl("out_any",fst_vals$File),"Out Any",
                                ifelse(grepl("out_across",fst_vals$File),
                                       "Out Combo", ifelse(grepl("within",fst_vals$File),"Out Within","Out All")))) 
fst_vals$Scenario<-ifelse(grepl("high",fst_vals$File),"High",
                           ifelse(grepl("extreme",fst_vals$File),"Extreme",
                                  ifelse(grepl("marg",fst_vals$File), "Marginal", 
                                         ifelse(grepl("panmictic",fst_vals$File),"Panmictic","Low")))) 
fst_vals$Linkage<-ifelse(grepl("lowlink",fst_vals$File),"Low Linkage","High Linkage")
fst_vals<-fst_vals[!grepl("4817126967349",fst_vals$File),]

View(fst_vals %>% filter(Scenario == "Panmictic") %>%
       filter(Linkage == "High Linkage") %>% pairwise_wilcox_test(
         V1 ~ Filter, paired = TRUE,
         p.adjust.method = "BH"
       ))


library(ggplot2)
library(cowplot)
fst_vals$Scenario = factor(fst_vals$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))
fst_vals %>% filter(Linkage == "Low Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
    ggplot(aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Fst") 

fst_vals %>% filter(Linkage == "High Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
    ggplot(aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Fst") 

csvs_fis<-csvs[grepl("*summary_stats_perloc*",csvs)]

fis<-lapply(csvs_fis,read.csv)
names(fis)<-csvs_fis
fis_vals<-do.call("rbind",fis)
fis_vals<-as.data.frame(fis_vals)
fis_vals$File<-rownames(fis_vals)
fis_vals<-fis_vals[!grepl("4817126967349",fis_vals$File),]
fis_vals$Filter<-ifelse(grepl("nohwe",fis_vals$File),"No Filter",
                         ifelse(grepl("out_any",fis_vals$File),"Out Any",
                                ifelse(grepl("out_across",fis_vals$File),
                                       "Out Combo", ifelse(grepl("within",fis_vals$File),"Out Within","Out All")))) 
fis_vals$Scenario<-ifelse(grepl("high",fis_vals$File),"High",
                           ifelse(grepl("extreme",fis_vals$File),"Extreme",
                                  ifelse(grepl("marg",fis_vals$File), "Marginal", 
                                         ifelse(grepl("panmictic",fis_vals$File),"Panmictic","Low")))) 
fis_vals$Linkage<-ifelse(grepl("lowlink",fis_vals$File),"Low Linkage","High Linkage")

fis_vals_sub<-fis_vals[fis_vals$Scenario=="Extreme",]

fis_vals_sub %>% filter(Linkage == "Low Linkage") %>% 
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
    ggplot(aes(x=Fis,y=Filter)) + geom_point() +
  facet_grid(.~Scenario,scales="free") + theme_classic() + ylab("HWE Filter") + xlab("Fis") + ylab("Fst") 

fis_cM<-lapply(fis,function(x){colMeans(x[2:ncol(x)],na.rm = T)})
fis_cM<-do.call("rbind",fis_cM)
fis_cM<-as.data.frame(fis_cM)
fis_cM$File<-rownames(fis_cM)
fis_cM<-fis_cM[!grepl("4817126967349",fis_cM$File),]
fis_cM$Filter<-ifelse(grepl("nohwe",fis_cM$File),"No Filter",
                      ifelse(grepl("out_any",fis_cM$File),"Out Any",
                             ifelse(grepl("out_across",fis_cM$File),
                                    "Out Combo", ifelse(grepl("within",fis_cM$File),"Out Within","Out All")))) 
fis_cM$Scenario<-ifelse(grepl("high",fis_cM$File),"High",
                        ifelse(grepl("extreme",fis_cM$File),"Extreme",
                               ifelse(grepl("marg",fis_cM$File), "Marginal", 
                                      ifelse(grepl("panmictic",fis_cM$File),"Panmictic","Low")))) 
fis_cM$Linkage<-ifelse(grepl("lowlink",fis_cM$File),"Low Linkage","High Linkage")

fis_cM %>% filter(Linkage == "High Linkage") %>% #filter(Filter=="HWE Out Combo") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
    ggplot(aes(x=Fst,y=Fis,colour=Filter, fill=Filter))+ geom_point() + #theme_minimal() +
  facet_wrap(Filter~Scenario,scales="free",ncol=5) + theme_bw() + ylab("Fis") + xlab("Fst") #+




fis_cM %>% filter(Linkage == "Low Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
    ggplot(aes(x=Hs,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Ho")


fis_regressions<-lapply(fis,function(x){t<-summary(lm(1-(x$Ho/x$Ht)~x$Fst))
t$coefficients[2,c(1,4)]
})


fis_regressions<-do.call("rbind",fis_regressions)
fis_regressions<-as.data.frame(fis_regressions)
colnames(fis_regressions)<-c("Slope","PVal")
fis_regressions$File<-rownames(fis_regressions)
fis_regressions<-fis_regressions[!grepl("4817126967349",fis_regressions$File),]
fis_regressions$Filter<-ifelse(grepl("nohwe",fis_regressions$File),"No Filter",
                               ifelse(grepl("out_any",fis_regressions$File),"Out Any",
                                      ifelse(grepl("out_across",fis_regressions$File),
                                             "Out Combo", ifelse(grepl("within",fis_regressions$File),"Out Within","Out All")))) 
fis_regressions$Scenario<-ifelse(grepl("high",fis_regressions$File),"High",
                                 ifelse(grepl("extreme",fis_regressions$File),"Extreme",
                                        ifelse(grepl("marg",fis_regressions$File), "Marginal", 
                                               ifelse(grepl("panmictic",fis_regressions$File),"Panmictic","Low")))) 
fis_regressions$Linkage<-ifelse(grepl("lowlink",fis_regressions$File),"Low Linkage","High Linkage")
fis_regressions$Scenario = factor(fis_regressions$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))

library(rstatix)            
View(fis_regressions %>% filter(Scenario == "Extreme") %>%
       filter(Linkage == "Low Linkage") %>% pairwise_wilcox_test(
  Slope ~ Filter, paired = TRUE,
  p.adjust.method = "BH"
))



#fis_regressions$PVal<-fis_regressions$PVal*10
#fis_regressions<-fis_regressions[fis_regressions$PVal<0.05,]





###### REAL DATA #########

### ISOPOD ####
library(stringr)
setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/real_data_within/isopod/")
isopod_PS_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
isopod_PS_names<-isopod_PS_names[grepl("PCA",isopod_PS_names)]
isopod_PS <- lapply(isopod_PS_names,read.csv)
names(isopod_PS) <- isopod_PS_names# paste0(word(c,1,2,sep="_"))
isopod_PS<-isopod_PS[grepl("*subrep*",names(isopod_PS))]
popmap<-read.csv("./isopod_popmap.txt")

file <- read.vcfR("./within_isopodhwe_out_within_subrep1.recode.vcf",verbose = FALSE) 
file_gl <- vcfR2genlight(file)
file_gl$pop<-as.factor(popmap$Pop[match(file_gl$ind.names,popmap$Ind)])
pops <- as.character(unique(popmap$Pop))
isopod_pops<-file_gl$pop
isopod_PS_within<-isopod_PS[grepl("within",names(isopod_PS))]
pcst_isopodPS_within <- lapply(isopod_PS_within,pc_st,isopod_pops)

file <- read.vcfR("./nofilt_isopod_genfilt_subrep1.recode.vcf",verbose = FALSE) 
file_gl <- vcfR2genlight(file)
file_gl$pop<-as.factor(popmap$Pop[match(file_gl$ind.names,popmap$Ind)])
pops <- as.character(unique(popmap$Pop))
isopod_pops<-file_gl$pop
isopod_PS<-isopod_PS[!grepl("within",names(isopod_PS))]
pcst_isopodPS <- lapply(isopod_PS,pc_st,isopod_pops)

pcst_isopodPS<-c(pcst_isopodPS,pcst_isopodPS_within)

pcst_isopodPS<-as.data.frame(do.call("rbind",pcst_isopodPS))
pcst_isopodPS$Filter<-rownames(pcst_isopodPS)

pcst_isopodPS$Scenario<-ifelse(grepl("within",pcst_isopodPS$Filter),"Out Within",
                           ifelse(grepl("all",pcst_isopodPS$Filter),"Out All",
                                  ifelse(grepl("any",pcst_isopodPS$Filter), "Out Any", 
                                         ifelse(grepl("across",pcst_isopodPS$Filter),"Out Combo","No Filter")))) 
isopod_Fst_names <- list.files(path = "./",pattern="*.csv",
                               recursive = T,full.names = T)
isopod_Fst_names<-isopod_Fst_names[grepl("Fst",isopod_Fst_names)]
isopod_Fst <- lapply(isopod_Fst_names,read.csv)
names(isopod_Fst) <- isopod_Fst_names# paste0(word(c,1,2,sep="_"))
isopod_Fst<-isopod_Fst[grepl("*subrep*",names(isopod_Fst))]

fst_vals_isopod <- lapply(isopod_Fst,function(x){mean(as.matrix(x[2:ncol(x)]),na.rm=T)})
fst_vals_isopod<-do.call("rbind",fst_vals_isopod)
fst_vals_isopod<-as.data.frame(fst_vals_isopod)
fst_vals_isopod$File<-rownames(fst_vals_isopod)
fst_vals_isopod$Scenario<-ifelse(grepl("within",fst_vals_isopod$File),"Out Within",
                                 ifelse(grepl("all",fst_vals_isopod$File),"Out All",
                                        ifelse(grepl("any",fst_vals_isopod$File), "Out Any", 
                                               ifelse(grepl("across",fst_vals_isopod$File),"Out Combo","No Filter")))) 

### ZEBRA ####
library(stringr)
setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/real_data_within/zebra/")
zebra_PS_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
zebra_PS_names<-zebra_PS_names[grepl("PCA",zebra_PS_names)]
zebra_PS <- lapply(zebra_PS_names,read.csv)
names(zebra_PS) <- zebra_PS_names# paste0(word(c,1,2,sep="_"))
zebra_PS<-zebra_PS[grepl("*subrep*",names(zebra_PS))]


pop_list<-read.csv("./zebmetdat.csv",head=T)
colnames(pop_list)<-c("V1","V2")


file <- read.vcfR("./within_zebrahwe_out_within_subrep1.recode.vcf",verbose = FALSE) 
file_gl <- vcfR2genlight(file)
file_gl$pop<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])
pops <- as.character(unique(pop_list$V2))
zebra_pops<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])
zebra_PS_within<-zebra_PS[grepl("within",names(zebra_PS))]
pcst_zebraPS_within <- lapply(zebra_PS_within,pc_st,zebra_pops)

file <- read.vcfR("./nofilt_populations.snps_subrep3.recode.vcf",verbose = FALSE) 
file_gl <- vcfR2genlight(file)
file_gl$pop<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])
pops <- as.character(unique(pop_list$V2))

zebra_pops<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])


zebra_PS<-zebra_PS[!grepl("within",names(zebra_PS))]
pcst_zebraPS <- lapply(zebra_PS,pc_st,zebra_pops)


pcst_zebraPS<-c(pcst_zebraPS,pcst_zebraPS_within)
names_zeb<-names(pcst_zebraPS)
pcst_zebraPS<-as.data.frame(do.call("rbind",pcst_zebraPS))
pcst_zebraPS$Filter<-names_zeb


pcst_zebraPS$Scenario<-ifelse(grepl("within",pcst_zebraPS$Filter),"Out Within",
                               ifelse(grepl("all",pcst_zebraPS$Filter),"Out All",
                                      ifelse(grepl("any",pcst_zebraPS$Filter), "Out Any", 
                                             ifelse(grepl("across",pcst_zebraPS$Filter),"Out Combo","No Filter")))) 

zebra_Fst_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
zebra_Fst_names<-zebra_Fst_names[grepl("Fst",zebra_Fst_names)]
zebra_Fst <- lapply(zebra_Fst_names,read.csv)
names(zebra_Fst) <- zebra_Fst_names# paste0(word(c,1,2,sep="_"))
zebra_Fst<-zebra_Fst[grepl("*subrep*",names(zebra_Fst))]

fst_vals_zebra <- lapply(zebra_Fst,function(x){mean(as.matrix(x[2:ncol(x)]),na.rm=T)})
fst_vals_zebra<-do.call("rbind",fst_vals_zebra)
fst_vals_zebra<-as.data.frame(fst_vals_zebra)
fst_vals_zebra$File<-rownames(fst_vals_zebra)
fst_vals_zebra$Scenario<-ifelse(grepl("within",fst_vals_zebra$File),"Out Within",
                                ifelse(grepl("all",fst_vals_zebra$File),"Out All",
                                       ifelse(grepl("any",fst_vals_zebra$File), "Out Any", 
                                              ifelse(grepl("across",fst_vals_zebra$File),"Out Combo","No Filter")))) 

### FURSEAL ####
library(stringr)
setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/real_data_within/furseal/")
seal_PS_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
seal_PS_names<-seal_PS_names[grepl("PCA",seal_PS_names)]
seal_PS <- lapply(seal_PS_names,read.csv)
names(seal_PS) <- seal_PS_names# paste0(word(c,1,2,sep="_"))
seal_PS<-seal_PS[grepl("*subrep*",names(seal_PS))]

file <- read.vcfR("./within_fursealhwe_out_within_subrep1.recode.vcf",verbose = FALSE) 
file_gl <- vcfR2genlight(file)
pop_list<-read.table("./popmap.tsv")
file_gl$pop<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])
pops <- as.character(unique(pop_list$V2))
seal_pops<-file_gl$pop
seal_PS_within<-seal_PS[grepl("within",names(seal_PS))]
pcst_sealPS_within <- lapply(seal_PS_within,pc_st,seal_pops)

file <- read.vcfR("./hwe_out_across_populations.snps_subrep9.recode.vcf",verbose = FALSE) 
file_gl <- vcfR2genlight(file)
file_gl$pop<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])
pops <- as.character(unique(pop_list$V2))
seal_pops<-file_gl$pop
seal_PS<-seal_PS[!grepl("within",names(seal_PS))]
pcst_sealPS<- lapply(seal_PS,pc_st,seal_pops)



pcst_sealPS<-c(pcst_sealPS,pcst_sealPS_within)

pcst_sealPS<-as.data.frame(do.call("rbind",pcst_sealPS))
pcst_sealPS$Filter<-rownames(pcst_sealPS)

pcst_sealPS$Scenario<-ifelse(grepl("within",pcst_sealPS$Filter),"Out Within",
                               ifelse(grepl("all",pcst_sealPS$Filter),"Out All",
                                      ifelse(grepl("any",pcst_sealPS$Filter), "Out Any", 
                                             ifelse(grepl("across",pcst_sealPS$Filter),"Out Combo","No Filter")))) 
library(stringr)
seal_Fst_names <- list.files(path = "./",pattern="*.csv",
                             recursive = T,full.names = T)
seal_Fst_names<-seal_Fst_names[grepl("Fst",seal_Fst_names)]
seal_Fst <- lapply(seal_Fst_names,read.csv)
names(seal_Fst) <- seal_Fst_names# paste0(word(c,1,2,sep="_"))
seal_Fst<-seal_Fst[grepl("*subrep*",names(seal_Fst))]

fst_vals_seal <- lapply(seal_Fst,function(x){mean(as.matrix(x[2:ncol(x)]),na.rm=T)})
fst_vals_seal<-do.call("rbind",fst_vals_seal)
fst_vals_seal<-as.data.frame(fst_vals_seal)
fst_vals_seal$File<-rownames(fst_vals_seal)
fst_vals_seal$Scenario<-ifelse(grepl("within",fst_vals_seal$File),"Out Within",
                               ifelse(grepl("all",fst_vals_seal$File),"Out All",
                                      ifelse(grepl("any",fst_vals_seal$File), "Out Any", 
                                             ifelse(grepl("across",fst_vals_seal$File),"Out Combo","No Filter")))) 

### JOINT ANALYSIS ####
real_pcst<-rbind(pcst_sealPS,pcst_isopodPS,pcst_zebraPS)
real_pcst$DatType<-c(rep("seal",50),rep("isopod",50),rep("zebra",50))
real_pcst$Scen_F = factor(real_pcst$DatType, levels=c('seal','zebra','isopod'))
(real_pcst %>% filter(Scen_F == "seal") %>% pairwise_wilcox_test(
  V1 ~ Scenario, paired = TRUE, exact=F,
  p.adjust.method = "BH"
))








real_fst<-rbind(fst_vals_seal,fst_vals_isopod,fst_vals_zebra)
real_fst$DatType<-c(rep("seal",50),rep("isopod",50),rep("zebra",50))
real_fst$Scen_F = factor(real_fst$DatType, levels=c('seal','zebra','isopod'))
real_fst$Scenario[real_fst$Scenario=="Out Across"] <- "Out Combo"
View(real_fst %>% filter(Scen_F == "seal") %>% pairwise_wilcox_test(
  V1 ~ Scenario, paired = TRUE, exact=F,
  p.adjust.method = "BH"
))




###### STRUCTURE
setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/structure_dat/")
k6_files<-list.files(pattern="*out_f")
k6_files<-k6_files[grepl("K6",k6_files)]
k6_files<-k6_files[!grepl("4817126967349",k6_files)]
k6_dat<-lapply(k6_files,read_struc_nucmat)
k6_dat<-do.call("rbind",k6_dat)
k6_dat<-data.frame(k6_dat,k6_files)
k6_dat$filter<-ifelse(grepl("nohwe",k6_dat$k6_files),"No Filter",
                          ifelse(grepl("out_any",k6_dat$k6_files),"Out Any",
                                 ifelse(grepl("out_across",k6_dat$k6_files),
                                        "Out Combo", ifelse(grepl("within",k6_dat$k6_files), 
                                                                 "Out Within", "Out All"))))
k6_dat$Scenario<-ifelse(grepl("high",k6_dat$k6_files),"High",
                                 ifelse(grepl("extreme",k6_dat$k6_files),"Extreme",
                                        ifelse(grepl("marg",k6_dat$k6_files), "Marginal", 
                                               ifelse(grepl("panmictic",k6_dat$k6_files),"Panmictic","Low")))) 
k6_dat$Linkage<-ifelse(grepl("lowlink",k6_dat$k6_files),"Low Linkage","High Linkage")
k6_dat$Scenario = factor(k6_dat$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))                      
k6_dat$Seed<-gsub('.*_([0-9]+).*','\\1',k6_dat$k6_files)
k6_dat<-k6_dat[!grepl(".gz",k6_dat$k6_files),]
struc_dat<-k6_dat
library(rstatix)
View(struc_dat %>% filter(Scenario == "Extreme") %>%
       filter(Linkage == "High Linkage") %>% pairwise_wilcox_test(
          k6_dat ~ filter, paired = TRUE,
         p.adjust.method = "BH"
       ))

library(tidyverse);library(ggridges)

neis_dist<-function(file){
  a<-read.csv(file)
  a$x
}
setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/pca_structure_output/")
neis_files<-list.files(pattern="*neis_dist.csv",recursive = T,full.names = T)
neis_files<-neis_files[grepl("subrep",neis_files)]
neis_files<-neis_files[!grepl("4817126967349",neis_files)]
neis_dat<-lapply(neis_files,neis_dist)
names(neis_dat)<-neis_files
neis_dat<-do.call("rbind",neis_dat)
neis_dat<-data.frame(rownames(neis_dat),neis_dat)
colnames(neis_dat)<-c("File","Neis")
neis_dat$filter<-ifelse(grepl("nohwe",neis_dat$File),"No Filter",
                      ifelse(grepl("out_any",neis_dat$File),"Out Any",
                             ifelse(grepl("out_across",neis_dat$File),
                                    "Out Combo", ifelse(grepl("within",neis_dat$File), 
                                                        "Out Within", "Out All"))))
neis_dat$Scenario<-ifelse(grepl("high",neis_dat$File),"High",
                        ifelse(grepl("extreme",neis_dat$File),"Extreme",
                               ifelse(grepl("marg",neis_dat$File), "Marginal", 
                                      ifelse(grepl("panmictic",neis_dat$File),"Panmictic","Low")))) 
neis_dat$Linkage<-ifelse(grepl("lowlink",neis_dat$File),"Low Linkage","High Linkage")
neis_dat$Scenario = factor(neis_dat$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))                      
neis_dat$Seed<-gsub('.*_([0-9]+).*','\\1',neis_dat$File)
neis_dat$Subrep<-word(neis_dat$File,1,sep="neis")
neis_dat$Subrep<-word(neis_dat$Subrep,-1,sep="_")
neis_dat$Shared<-paste(neis_dat$Seed,neis_dat$filter,neis_dat$Scenario,neis_dat$Linkage,neis_dat$Subrep)
k6_dat$Subrep<-word(k6_dat$k6_files,2,sep="[.]")
k6_dat$Subrep<-word(k6_dat$Subrep,-1,sep="_")
k6_dat$Shared<-paste(k6_dat$Seed,k6_dat$filter,k6_dat$Scenario,k6_dat$Linkage,k6_dat$Subrep)
k6_dat<-k6_dat[!grepl("*gz",k6_dat$k6_files),]
neis_dat$Struc<-k6_dat$k6_dat[match(neis_dat$Shared,k6_dat$Shared)]
neis_dat$Agg_Var<-paste(neis_dat$Seed,neis_dat$filter,neis_dat$Scenario,neis_dat$Linkage)
coeffs_neis_struc<-lmList(Struc ~ Neis | Agg_Var, data=neis_dat)
coefs<-lapply(coeffs_neis_struc,function(x){x$coefficients})
coefs<-do.call("rbind",coefs)
coefs<-data.frame(coefs)
temp<-names(coeffs_neis_struc)
temp<-as.data.frame(do.call("rbind",str_split(temp,pattern = " ")))
temp$Filter<-paste(temp$V2,temp$V3);temp$Linkage<-paste(temp$V5,temp$V6)
temp$V2<-NULL;temp$V3<-NULL;temp$V5<-NULL;temp$V6<-NULL
coefs<-cbind(coefs,temp)
colnames(coefs)[3:4]<-c("Seed","PopStruc")

coefs %>% filter(Linkage == "Low Linkage") %>%
  mutate(filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=X.Intercept.,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~PopStruc,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Slow")


setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/structure_dat/")
k6_files<-list.files(pattern="*out_f")
k6_files<-k6_files[grepl("K6",k6_files)]
k6_files<-k6_files[!grepl("4817126967349",k6_files)]
library(tidyverse)
variance_ln<-function(file){
  z<-readLines(file)
  z<-z[grepl("Mean value of ln likelihood",z)]
  as.numeric(word(z,2,sep="= "))
}
variances<-lapply(k6_files,variance_ln)
names(variances)<-k6_files
variances<-do.call("rbind",variances)
variances<-data.frame(variances,k6_files)
variances$filter<-ifelse(grepl("nohwe",variances$k6_files),"No Filter",
                      ifelse(grepl("out_any",variances$k6_files),"Out Any",
                             ifelse(grepl("out_across",variances$k6_files),
                                    "Out Combo", ifelse(grepl("within",variances$k6_files), 
                                                        "Out Within", "Out All"))))
variances$Scenario<-ifelse(grepl("high",variances$k6_files),"High",
                        ifelse(grepl("extreme",variances$k6_files),"Extreme",
                               ifelse(grepl("marg",variances$k6_files), "Marginal", 
                                      ifelse(grepl("panmictic",variances$k6_files),"Panmictic","Low")))) 
variances$Linkage<-ifelse(grepl("lowlink",variances$k6_files),"Low Linkage","High Linkage")
variances$Scenario = factor(variances$Scenario, levels=c('Panmictic','Marginal','Low','High','Extreme'))                      
variances$Seed<-gsub('.*_([0-9]+).*','\\1',variances$k6_files)
variances<-variances[!grepl("*gz",variances$k6_files),]
variances$VarStandVariable<-paste(variances$Seed,variances$Linkage,variances$filter,variances$Scenario,sep="_")
stats::aggregate(.~variances$VarStandVariable,variances$variances,mean)
maxed<-aggregate(variances$variances, list(variances$VarStandVariable), max)
mined<-aggregate(variances$variances, list(variances$VarStandVariable), min)
maxed$Dif<-maxed$x-mined$x
library(ggridges)
variances$PerGroup_MaxVar<-maxed$Dif[match(variances$VarStandVariable,maxed$Group.1)]
variances$StandVar<-variances$variances/variances$PerGroup_MaxVar
struc_variance_plot<- variances %>% filter(Linkage == "High Linkage") %>%
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=variances,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Average Nucleotide Distance")
struc_variance_plot

#ggsave(struc_variance_plot,filename = "/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/Figures/variance_plots_LowLinkage.pdf",width = 11.6,height=3.86)
#ggsave(struc_variance_plot,filename = "/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/Figures/variance_plots_HighLinkage.pdf",width = 11.6,height=3.86)

## REAL STRUC VARIANCE
setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/")
files<-list.files(pattern="*out_f",recursive = T)
files<-files[grepl("isopod|seal|zebra",files)]
files<-files[!grepl("structure_dat",files)]

library(tidyverse)
variance_ln<-function(file){
  z<-readLines(file)
  z<-z[grepl("Mean value of ln likelihood",z)]
  as.numeric(word(z,2,sep="= "))
}
variances<-lapply(files,variance_ln)
names(variances)<-files
variances<-do.call("rbind",variances)
variances<-data.frame(variances,files)
variances$filter<-ifelse(grepl("nofilt",variances$files),"No Filter",
                         ifelse(grepl("out_any",variances$files),"Out Any",
                                ifelse(grepl("out_across",variances$files),
                                       "Out Combo", ifelse(grepl("within",variances$files), 
                                                           "Out Within", "Out All"))))

variances$Scen_F = c(rep("Isopod",50),rep("Seal",50),rep("Zebra",50))
library(ggridges)
struc_variance_plot<- variances %>% 
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  mutate(Scen_F = fct_relevel(Scen_F, levels = "Seal","Zebra","Isopod")) %>%
  ggplot(aes(x=variances,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Average Nucleotide Distance")
struc_variance_plot

ggsave(struc_variance_plot,filename = "/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/Figures/variances_realdat.pdf",width = 11.6,height=3.86)

######### REAL STRUC ########

setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/seal_structure/")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k8)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                      ifelse(grepl("out_any",struc_mats$files),"Out Any",
                             ifelse(grepl("out_across",struc_mats$files),
                                    "Out Combo", ifelse(grepl("within",struc_mats$files), 
                                                        "Out Within", "Out All"))))

struc_mats_seal <-struc_mats
library(dplyr)
struc_mats_seal%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))

setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/zebra_structure//")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k9)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "Out Combo", ifelse(grepl("within",struc_mats$files), 
                                                            "Out Within", "Out All"))))
struc_mats_zebra<-struc_mats
library(dplyr)
struc_mats_zebra%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))


setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/isopod_structure//")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k8)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "Out Combo", ifelse(grepl("within",struc_mats$files), 
                                                            "Out Within", "Out All"))))
library(dplyr)
struc_mats%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))

struc_mats_isopod <- struc_mats
struc_real<- do.call("rbind",list(struc_mats_seal,struc_mats_zebra,struc_mats_isopod))
struc_real$Scen_F = c(rep("Seal",50),rep("Zebra",50),rep("Isopod",50))


(struc_real %>% filter(Scen_F == "Seal") %>% pairwise_wilcox_test(
  struc_mats ~ filter, paired = TRUE, exact=F,
  p.adjust.method = "BH"
))


######## Liklihoo