library(dartR)
library(vcfR)
library(devtools)
library(stringr)
#library(radiator)
setwd("/home/peawi142/HWE_Simulations/real_data/seals/single_SNPs/")
seal_PS_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
seal_PS_names<-seal_PS_names[grepl("PCA",seal_PS_names)]
seal_PS <- lapply(seal_PS_names,read.csv)
a<-paste(word(seal_PS_names,3,sep="/"))
b<-paste(word(seal_PS_names,3,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(seal_PS) <- paste0(a)
seal_PS<-seal_PS[grepl("*subrep*",names(seal_PS))]
seal_pops<-read.table("./popmap.tsv")
seal_pops<-seal_pops$V2
pcst_sealPS <- lapply(seal_PS,pc_st,seal_pops)
pcst_sealPS<-as.data.frame(do.call("rbind",pcst_sealPS))
pcst_sealPS$Filter<-rownames(pcst_sealPS)
pcst_sealPS$Filter <- word(pcst_sealPS$Filter,3,4,sep="_")
pcst_sealPS$Filter <- ifelse(grepl("*nofilt*",rownames(pcst_sealPS)),"nofilt",
                               word(rownames(pcst_sealPS),1,sep="_"))
ggplot(pcst_sealPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_sealPS$Filter<-ifelse(pcst_sealPS$Filter=="nofilt","No Filt",
                             ifelse(pcst_sealPS$Filter=="any","HWE Out Any",
                                    ifelse(pcst_sealPS$Filter=="all","HWE Out All","HWE Out Across")))


library(ggridges)
ggplot(pcst_sealPS,aes(x=V1,y=Filter)) + theme_minimal() + ggtitle("specified bandwidth of 0.01") +
  xlab("PCst") + stat_density_ridges(bandwidth = 0.01,stat="density")
ggplot(pcst_sealPS[pcst_sealPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("seal_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_sealPS$V1[pcst_sealPS$Filter=="HWE Out All"],
       pcst_sealPS$V1[pcst_sealPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_sealPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()


fst_simulations<-list.files(pattern="*csv",recursive = T,full.names = T)
fst_simulations <- fst_simulations[grepl("Fst",fst_simulations)]
fst_simulations_dat<-lapply(fst_simulations,read.csv)
names(fst_simulations_dat)<-fst_simulations
fst_simulations_dat<-fst_simulations_dat[grepl("subrep",names(fst_simulations_dat))]
a<-vector()
b<-vector()
c<-vector()
for ( i in 1:length(fst_simulations_dat)){
  t<-names(fst_simulations_dat)
  t<-as.character(word(t[i],1,sep="_"))
  #a[i]<-(mean(as.numeric(as.matrix(fst_filesdat[[t]])),na.rm=T))
  b[i]<-mean(as.matrix(fst_simulations_dat[[i]][2:9]),na.rm=T)
  c[i]<-t
}
seal_fsts<-data.frame(b,c)
colnames(seal_fsts)<-c("Fst","Filter")
seal_fsts$Filter<-ifelse(seal_fsts$Filter=="./nofilt","No Filt",
                           ifelse(seal_fsts$Filter=="./any","HWE Out Any",
                                  ifelse(seal_fsts$Filter=="./all","HWE Out All","HWE Out Across")))
ggplot(seal_fsts,aes(x=Fst,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("seal_PS_Fst") +
  xlab("Fst")


# summary stats

setwd("/home/peawi142/HWE_Simulations/real_data/seals/")
files_stat_names <- list.files(path = "./",pattern="*summary_stats.csv",
                               recursive = T,full.names = T)
files_stat_names<-files_stat_names[grepl("subrep",files_stat_names)]
files_stat<-lapply(files_stat_names,read.csv)
extract_fis<-function(file){
  file$x
}
seal_fis_dat<-lapply(files_stat,extract_fis)
newnames<-word(files_stat_names,4,sep="/") %>%
  word(.,1,2,sep="_") 
seal_fis_dat<-data.frame(do.call("rbind",seal_fis_dat));colnames(seal_fis_dat)<-c("Ho","Hs","Ht","Dst","Dstp","Fst","Fstp","Fis","Dest")
seal_fis_dat$newname<-newnames
seal_fis_dat$Filt<-word(newnames,1,sep="_")
seal_fis_distr <- ggplot(seal_fis_dat,aes(x=Ho,y=Filt)) + theme_minimal() + xlab("Ho") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1) + ggtitle("seal")
seal_fis_distr
####### STRUCTURE
setwd("~/HWE_Simulations/real_data/seals/structure/")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k8)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"HWE Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "HWE Out Across", "HWE Out All"))) 
struc_mats_seal<-struc_mats
library(dplyr)
struc_mats_seal%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))


ggplot(struc_mats_seal,aes(x=struc_mats,y=filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
  ) +
  xlab("Structure Distance")#+xlim(-0.002,0.002)
