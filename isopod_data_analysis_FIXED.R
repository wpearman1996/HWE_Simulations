library(dartR)
library(vcfR)
library(devtools)
#library(radiator)
setwd("/home/peawi142/HWE_Simulations/real_data/isopods/vcfs/VCFS_r2/monomorph_fix/ploidy_HWE_nloci_testing/poly_vcfs/")

isopod_PS_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
isopod_PS_names<-isopod_PS_names[grepl("PCA",isopod_PS_names)]
isopod_PS <- lapply(isopod_PS_names,read.csv)
a<-paste(word(isopod_PS_names,3,sep="/"))
b<-paste(word(isopod_PS_names,3,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(isopod_PS) <- paste0(word(c,1,2,sep="_"))
isopod_PS<-isopod_PS[grepl("*subrep*",names(isopod_PS))]
file <- read.vcfR("../../../nofilt_fixed_isopod_vcf_subrep5.recode.vcf",verbose = FALSE)
file_gl <- vcfR2genlight(file)

popmap<-read.csv("~/HWE_Simulations/real_data/isopods/metadata.csv")
file_gl$pop<-as.factor(popmap$pop[match(file_gl$ind.names,popmap$id)])
file_gl$pop[64]<-"Mahia Peninsula"
file_gl$pop[79]<-"Worser Bay"

file_gl$pop[111]<-"Mt Maunganui"

file_gl$pop[114]<-"Mt Maunganui"
pop<-as.character(file_gl$pop)
file_gl$pop<-as.factor(as.character((ifelse(pop=="Stanmore Bay Old","Auckland",
                                            ifelse(pop=="StanmoreBay", "Auckland",
                                                   ifelse(pop=="Hatfields Beach","Auckland",
                                                          ifelse(pop=="Browns Bay","Auckland",pop)))))))
isopod_pops<-file_gl$pop
pcst_isopodPS <- lapply(isopod_PS,pc_st,isopod_pops)
pcst_isopodPS<-as.data.frame(do.call("rbind",pcst_isopodPS))
pcst_isopodPS$Filter<-rownames(pcst_isopodPS)
#pcst_isopodPS$Filter <- word(pcst_isopodPS$Filter,3,4,sep="_")
pcst_isopodPS$Filter <- word(rownames(pcst_isopodPS),1,sep="_")
ggplot(pcst_isopodPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_isopodPS$Filter<-ifelse(pcst_isopodPS$Filter=="nofilt","No Filt",
                             ifelse(pcst_isopodPS$Filter=="any","HWE Out Any",
                                    ifelse(pcst_isopodPS$Filter=="all","HWE Out All","HWE Out Across")))



ggplot(pcst_isopodPS[pcst_isopodPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("isopod_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_isopodPS$V1[pcst_isopodPS$Filter=="HWE Out All"],
       pcst_isopodPS$V1[pcst_isopodPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_isopodPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()
ggplot(pcst_isopodPS,aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
  ) + xlab("Inferred PCst")#+xlim(-0.002,0.002)


fst_simulations<-list.files(pattern="*csv",recursive = T,full.names = T)
fst_simulations <- fst_simulations[grepl("Fst",fst_simulations)]
fst_simulations_dat<-lapply(fst_simulations,read.csv)
names(fst_simulations_dat)<-fst_simulations
fst_simulations_dat<-fst_simulations_dat[grepl("subrep",names(fst_simulations_dat))]
a<-vector()
b<-vector()
c<-vector()
names(fst_simulations_dat)<-word(names(fst_simulations_dat),2,sep="/")
for ( i in 1:length(fst_simulations_dat)){
  t<-names(fst_simulations_dat)
  t<-as.character(paste0("sim_",word(t[i],1,2,sep="_")))
  #a[i]<-(mean(as.numeric(as.matrix(fst_filesdat[[t]])),na.rm=T))
  b[i]<-mean(as.matrix(fst_simulations_dat[[i]][2:9]),na.rm=T)
  c[i]<-t
}
isopod_fsts<-data.frame(b,c)
colnames(isopod_fsts)<-c("Fst","Filter")
isopod_fsts$Filter<-word(isopod_fsts$Filter,1,2,sep="_")
isopod_fsts$Filter<-ifelse(isopod_fsts$Filter=="sim_nofilt","No Filt",
                           ifelse(isopod_fsts$Filter=="sim_any","HWE Out Any",
                                  ifelse(isopod_fsts$Filter=="sim_all","HWE Out All","HWE Out Across")))
ggplot(isopod_fsts,aes(x=Fst,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("isopod_PS_Fst") +
  xlab("Fst")

ggplot(isopod_fsts,aes(x=Fst,y=Filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
                      bandwidth=0.006
  ) 


####### STRUCTURE
setwd("~/HWE_Simulations/real_data/isopods/ploidy_HWE_nloci_testing/structure/")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k8)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"HWE Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "HWE Out Across", "HWE Out All"))) 

library(dplyr)
struc_mats%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))


#struc_mats<-rbind(struc_mats,
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",],
#                  struc_mats[struc_mats$filter=="HWE Out Across",]) # we need to add pseudo values to make ridgeplot add a distibution - we will modify this for the figure later.

isopod_struc_nucdist <- ggplot(struc_mats,aes(x=struc_mats,y=filter)) + theme_minimal() + xlab("Nuc Dist") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1,bandwidth=0.01)
isopod_struc_mats <- struc_mats
isopod_struc_nucdist
### isopod with polyploidy
setwd("/home/peawi142/HWE_Simulations/real_data/isopods/vcfs/")

isopod_poly_PS_names <- list.files(path = "./",pattern="*.csv",
                                   recursive = F,full.names = T)
isopod_poly_PS_names<-isopod_poly_PS_names[grepl("PCA",isopod_poly_PS_names)]
isopod_poly_PS <- lapply(isopod_poly_PS_names,read.csv)
a<-paste(word(isopod_poly_PS_names,3,sep="/"))
b<-paste(word(isopod_poly_PS_names,3,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(isopod_poly_PS) <-c# paste0(word(c,1,2,sep="_"))
isopod_poly_PS<-isopod_poly_PS[grepl("*subrep*",names(isopod_poly_PS))]
popmap<-read.csv("~/HWE_Simulations/real_data/isopods/metadata.csv")
isopod_Pops_poly<-as.factor(popmap$pop[match(isopod_poly_PS$isopods_hwe_out_across_isopod_macfilt.recode_subrep1$X,popmap$id)])
pcst_isopod_polyPS <- lapply(isopod_poly_PS,pc_st,isopod_Pops_poly)
pcst_isopod_polyPS<-as.data.frame(do.call("rbind",pcst_isopod_polyPS))
pcst_isopod_polyPS$Filter<-rownames(pcst_isopod_polyPS)
#pcst_isopod_polyPS$Filter <- word(pcst_isopod_polyPS$Filter,3,4,sep="_")
pcst_isopod_polyPS$Filter <- word(rownames(pcst_isopod_polyPS),2,4,sep="_")
ggplot(pcst_isopod_polyPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_isopod_polyPS$Filter<-ifelse(pcst_isopod_polyPS$Filter=="hwe_out_across","Out Across",
                                  ifelse(pcst_isopod_polyPS$Filter=="hwe_out_any","Out Any",
                                         ifelse(pcst_isopod_polyPS$Filter=="hwe_out_all","Out All","No Filter")))


ggplot(pcst_isopod_polyPS,aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
  ) + xlab("PCst")#+xlim(-0.002,0.002)


fst_simulations_poly<-list.files(pattern="*csv",recursive = F,full.names = T)
fst_simulations_poly <- fst_simulations_poly[grepl("Fst",fst_simulations_poly)]
fst_simulations_poly_dat<-lapply(fst_simulations_poly,read.csv)
names(fst_simulations_poly_dat)<-fst_simulations_poly
fst_simulations_poly_dat<-fst_simulations_poly_dat[grepl("subrep",names(fst_simulations_poly_dat))]
a<-vector()
b<-vector()
c<-vector()
names(fst_simulations_poly_dat)<-word(names(fst_simulations_poly_dat),2,sep="/")
for ( i in 1:length(fst_simulations_poly_dat)){
  t<-names(fst_simulations_poly_dat)
  t<-as.character(paste0("sim_",word(t[i],3,4,sep="_")))
  #a[i]<-(mean(as.numeric(as.matrix(fst_filesdat[[t]])),na.rm=T))
  b[i]<-mean(as.matrix(fst_simulations_poly_dat[[i]][2:9]),na.rm=T)
  c[i]<-t
}
isopod_fsts_poly<-data.frame(b,c)
colnames(isopod_fsts_poly)<-c("Fst","Filter")
isopod_fsts_poly$Filter<-word(isopod_fsts_poly$Filter,2,3,sep="_")
isopod_fsts_poly$Filter<-ifelse(isopod_fsts_poly$Filter=="out_across","Out Across",
                                ifelse(isopod_fsts_poly$Filter=="out_any","Out Any",
                                       ifelse(isopod_fsts_poly$Filter=="out_all","Out All","No Filter")))

ggplot(isopod_fsts_poly,aes(x=Fst,y=Filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
                      bandwidth=0.006
  ) 


setwd("~/HWE_Simulations/real_data/isopods/structure/")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k8)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"HWE Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "HWE Out Across", "HWE Out All"))) 

library(dplyr)
struc_mats%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))


struc_mats<-rbind(struc_mats,
                  struc_mats[struc_mats$filter=="HWE Out Across",],
                  struc_mats[struc_mats$filter=="HWE Out Across",],
                  struc_mats[struc_mats$filter=="HWE Out Across",],
                  struc_mats[struc_mats$filter=="HWE Out Across",],
                  struc_mats[struc_mats$filter=="HWE Out Across",],
                  struc_mats[struc_mats$filter=="HWE Out Across",],
                  struc_mats[struc_mats$filter=="HWE Out Any",],
                  struc_mats[struc_mats$filter=="HWE Out Any",],
                  struc_mats[struc_mats$filter=="HWE Out Any",]
) # we need to add pseudo values to make ridgeplot add a distibution - we will modify this for the figure later.
struc_mats$filter<-gsub("HWE ","",struc_mats$filter)
isopod_polyploidy_struc_nucdist <- ggplot(struc_mats,aes(x=struc_mats,y=filter)) + theme_minimal() + xlab("Nuc Dist") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1,bandwidth=0.01)
isopod_polyploidy_struc_mats <- struc_mats

isopod_polyploidy_struc_mats_plot<-ggplot(isopod_polyploidy_struc_mats,aes(x=struc_mats,y=filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000000000001,
  ) +
  xlab("StructureNucDist")#+xlim(-0.002,0.002)

plot_grid(isopod_fst_means, isopod_pcst_fst,isopod_pcst_fst_noacross,
          isopod_pcst_all,isopod_pcst_noacross,isopod_struc_nucdist
          , labels = c('InfFst_Means isopod PS','PCst vs InfFst', 'PCst vs InfFst Zoomed','PCst',
                       'PCst zoomed', 'Structure Nuc Dist'))


# summary stats
setwd("/home/peawi142/HWE_Simulations/real_data/isopods/vcfs/VCFS_r2/")
files_stat_names <- list.files(path = "./",pattern="*summary_stats.csv",
                               recursive = F,full.names = T)
files_stat_names<-files_stat_names[grepl("subrep",files_stat_names)]
files_stat<-lapply(files_stat_names,read.csv)
extract_fis<-function(file){
  file$x
}
isopod_ploidy_fis_dat<-lapply(files_stat,extract_fis)
newnames<-word(files_stat_names,3,sep="/") %>%
  word(.,1,2,sep="_") 
isopod_ploidy_fis_dat<-data.frame(do.call("rbind",isopod_ploidy_fis_dat));colnames(isopod_ploidy_fis_dat)<-c("Ho","Hs","Ht","Dst","Dstp","Fst","Fstp","Fis","Dest")
isopod_ploidy_fis_dat$newname<-newnames
isopod_ploidy_fis_dat$Filt<-word(newnames,1,sep="_")

setwd("/home/peawi142/HWE_Simulations/real_data/isopods/vcfs/VCFS_r2/monomorph_fix/")
files_stat_names <- list.files(path = "./",pattern="*summary_stats.csv",
                               recursive = T,full.names = T)
files_stat_names<-files_stat_names[grepl("subrep",files_stat_names)]
files_stat<-lapply(files_stat_names,read.csv)
extract_fis<-function(file){
  file$x
}
isopod_fis_dat<-lapply(files_stat,extract_fis)
newnames<-word(files_stat_names,3,sep="/") %>%
  word(.,1,2,sep="_") 
isopod_fis_dat<-data.frame(do.call("rbind",isopod_fis_dat));colnames(isopod_fis_dat)<-c("Ho","Hs","Ht","Dst","Dstp","Fst","Fstp","Fis","Dest")
isopod_fis_dat$newname<-newnames
isopod_fis_dat$Filt<-word(newnames,1,sep="_")
isopod_fis_dat$DatType<-"NoPloidy"
isopod_ploidy_fis_dat$DatType<-"PolyPloid"
isopod_fis_dat <- rbind(isopod_fis_dat,isopod_ploidy_fis_dat)
ggplot(isopod_fis_dat,aes(x=Hs,y=Filt))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red", bandwidth=0.001) +
  theme_bw() + ylab("HWE Filter") + facet_grid(.~DatType,scales="free") 

