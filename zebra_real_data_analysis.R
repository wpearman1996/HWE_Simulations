library(dartR)
library(vcfR)
library(devtools)
library(stringr)
#library(radiator)
setwd("/home/peawi142/HWE_zebra/real_data/zebra/single_SNPs/")
zebra_PS_names <- list.files(path = "./",pattern="*.csv",
                              recursive = T,full.names = T)
zebra_PS_names<-zebra_PS_names[grepl("PCA",zebra_PS_names)]
zebra_PS <- lapply(zebra_PS_names,read.csv)
a<-paste(word(zebra_PS_names,3,sep="/"))
b<-paste(word(zebra_PS_names,3,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(zebra_PS) <- paste0(a)
zebra_PS<-zebra_PS[grepl("*subrep*",names(zebra_PS))]
zebra_pops<-read.table("./zebmetdat.csv")
zebra_pops<-zebra_pops$V2
file <- read.vcfR("nofilt_populations.snps_subrep10.recode.vcf",verbose = FALSE)
file_gl <- vcfR2genlight(file)

metadat<-read.csv("~/HWE_Simulations/real_data/zebra/single_SNPs/zebmetdat.csv",head=F)
colnames(metadat)<-c("Name","Locality")
pop_list<-metadat$Locality
zebra_pops<-as.factor(metadat$Locality[match(file_gl$ind.names,metadat$Name)])

pcst_zebraPS <- lapply(zebra_PS,pc_st,zebra_pops)
pcst_zebraPS<-as.data.frame(do.call("rbind",pcst_zebraPS))
pcst_zebraPS$Filter<-rownames(pcst_zebraPS)
pcst_zebraPS$Filter <- word(pcst_zebraPS$Filter,3,4,sep="_")
pcst_zebraPS$Filter <- ifelse(grepl("*nofilt*",rownames(pcst_zebraPS)),"nofilt",
                               word(rownames(pcst_zebraPS),1,sep="_"))
ggplot(pcst_zebraPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_zebraPS$Filter<-ifelse(pcst_zebraPS$Filter=="nofilt","No Filt",
                             ifelse(pcst_zebraPS$Filter=="any","HWE Out Any",
                                    ifelse(pcst_zebraPS$Filter=="all","HWE Out All","HWE Out Across")))


library(ggridges)
ggplot(pcst_zebraPS,aes(x=V1,y=Filter)) + theme_minimal() + ggtitle("specified bandwidth of 0.01") +
  xlab("PCst") + stat_density_ridges(bandwidth = 0.01,stat="density")
ggplot(pcst_zebraPS[pcst_zebraPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("zebra_PS_PCA_PCst") +
  xlab("PCst")

ggplot(pcst_zebraPS,aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
  ) + xlab("Inferred PCst")#+xlim(-0.002,0.002)

t.test(pcst_zebraPS$V1[pcst_zebraPS$Filter=="HWE Out All"],
       pcst_zebraPS$V1[pcst_zebraPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_zebraPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()


#setwd("./vcfs/")
fst_zebra<-list.files(pattern="*csv",recursive = T,full.names = T)
fst_zebra <- fst_zebra[grepl("Fst",fst_zebra)]
fst_zebra_dat<-lapply(fst_zebra,read.csv)
names(fst_zebra_dat)<-fst_zebra
fst_zebra_dat<-fst_zebra_dat[grepl("subrep",names(fst_zebra_dat))]
a<-vector()
b<-vector()
c<-vector()
for ( i in 1:length(fst_zebra_dat)){
  t<-names(fst_zebra_dat)
  t<-as.character(paste0(word(t[i],1,sep="_")))
  #a[i]<-(mean(as.numeric(as.matrix(fst_filesdat[[t]])),na.rm=T))
  b[i]<-mean(as.matrix(fst_zebra_dat[[i]][2:10]),na.rm=T)
  c[i]<-t
}
zebra_fsts<-data.frame(b,c)
colnames(zebra_fsts)<-c("Fst","Filter")
zebra_fsts$Filter<-ifelse(zebra_fsts$Filter=="./nofilt","No Filt",
                           ifelse(zebra_fsts$Filter=="./any","HWE Out Any",
                                  ifelse(zebra_fsts$Filter=="./all","HWE Out All","HWE Out Across")))
ggplot(zebra_fsts,aes(x=Fst,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("zebra_PS_Fst") +
  xlab("Fst")

ggplot(zebra_fsts, aes(x=Fst, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

# summary stats
setwd("/home/peawi142/HWE_Simulations/real_data/zebra/single_SNPs/")
files_stat_names <- list.files(path = "./",pattern="*summary_stats.csv",
                               recursive = T,full.names = T)
files_stat_names<-files_stat_names[grepl("subrep",files_stat_names)]
files_stat<-lapply(files_stat_names,read.csv)
extract_fis<-function(file){
  file$x
}
zebra_fis_dat<-lapply(files_stat,extract_fis)
newnames<-word(files_stat_names,3,sep="/") %>%
  word(.,1,2,sep="_") 
zebra_fis_dat<-data.frame(do.call("rbind",zebra_fis_dat));colnames(zebra_fis_dat)<-c("Ho","Hs","Ht","Dst","Dstp","Fst","Fstp","Fis","Dest")
zebra_fis_dat$newname<-newnames
zebra_fis_dat$Filt<-word(newnames,1,sep="_")
zebra_fis_distr <- ggplot(zebra_fis_dat,aes(x=Ho,y=Filt)) + theme_minimal() + xlab("Ho") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1) + ggtitle("zebra")
zebra_fis_distr


setwd("/home/peawi142/HWE_Simulations/real_data/zebras/")
files_stat_names <- list.files(path = "./",pattern="*summary_stats.csv",
                               recursive = T,full.names = T)
files_stat_names<-files_stat_names[grepl("subrep",files_stat_names)]
files_stat<-lapply(files_stat_names,read.csv)
extract_fis<-function(file){
  file$x
}
zebra_fis_dat<-lapply(files_stat,extract_fis)
newnames<-word(files_stat_names,4,sep="/") %>%
  word(.,1,2,sep="_") 
zebra_fis_dat<-data.frame(do.call("rbind",zebra_fis_dat));colnames(zebra_fis_dat)<-c("Ho","Hs","Ht","Dst","Dstp","Fst","Fstp","Fis","Dest")
zebra_fis_dat$newname<-newnames
zebra_fis_dat$Filt<-word(newnames,1,sep="_")
zebra_fis_distr <- ggplot(zebra_fis_dat,aes(x=Ho,y=Filt)) + theme_minimal() + xlab("Ho") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1) + ggtitle("zebra")
zebra_fis_distr
####### STRUCTURE
setwd("~/HWE_Simulations/real_data/zebra/structure/")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat_k9)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nofilt",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"HWE Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "HWE Out Across", "HWE Out All"))) 
struc_mats_zebra<-struc_mats
library(dplyr)
struc_mats_zebra%>%
  group_by(filter)%>% 
  summarise(Mean=mean(struc_mats), Max=max(struc_mats), Min=min(struc_mats), Median=median(struc_mats), Std=sd(struc_mats))


ggplot(struc_mats_zebra,aes(x=struc_mats,y=filter))+  #theme_minimal() +
  geom_density_ridges(quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      jittered_points = TRUE, alpha = 0.7,
                      vline_size = 1, vline_color = "red",
                      point_size = 0.4, point_alpha = 0.3,rel_min_height = 0.000001,
  ) +
  xlab("Structure Distance")#+xlim(-0.002,0.002)

