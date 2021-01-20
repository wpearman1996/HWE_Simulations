library(dartR)
library(vcfR)
library(devtools)
#library(radiator)
setwd("/home/peawi142/HWE_Simulations/")
setwd("/home/peawi142/HWE_Study_Dec2020/slim3.5_newsims/real_data/isopods/")
data <- gl.read.dart(filename = ".//Report_DIsoc19-4130_1_moreOrders_SNP_2.csv", 
                     ind.metafile = ".//metadata.csv", probar = TRUE)
iso_dat<-gl.drop.ind(data,ind.list = "BBFA3", recalc = TRUE)
iso_dat<-gl.filter.callrate(iso_dat,threshold = 0.8,method="loc")
iso_dat<-gl.filter.repavg(iso_dat,threshold=0.9)
iso_dat<-gl.filter.monomorphs(iso_dat)
iso_dat<-gl.filter.secondaries(iso_dat,method="best")
install_github("thierrygosselin/radiator")
library(radiator)
iso_rad<-tidy_genlight(data = iso_dat, tidy = TRUE, gds = TRUE)
iso_rad$GT_VCF<-ifelse(iso_rad$GT_BIN == 0,"0/0",
                       ifelse(iso_rad$GT_BIN == 1,"1/1",
                              ifelse(iso_rad$GT_BIN == 2,"0/1","./.")))
iso_rad$GT_VCF<-replace_na(iso_rad$GT_VCF,"./.")
genomic_converter(iso_rad,output = "vcf",filename="sdgs")

isopod_PS_names <- list.files(path = "./",pattern="*.csv",
                            recursive = T,full.names = T)
isopod_PS_names<-isopod_PS_names[grepl("PCA",isopod_PS_names)]
isopod_PS <- lapply(isopod_PS_names,read.csv)
a<-paste(word(isopod_PS_names,3,sep="/"))
b<-paste(word(isopod_PS_names,3,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(isopod_PS) <- paste0(a,c)
isopod_PS<-isopod_PS[grepl("*subrep*",names(isopod_PS))]
isopod_pops<-readLines("~/HWE_Study_Dec2020/slim3.5_newsims/real_data/isopods/isopod_popmap.txt")
pcst_isopodPS <- lapply(isopod_PS,pc_st,isopod_pops)
pcst_isopodPS<-as.data.frame(do.call("rbind",pcst_isopodPS))
pcst_isopodPS$Filter<-rownames(pcst_isopodPS)
pcst_isopodPS$Filter <- word(pcst_isopodPS$Filter,3,4,sep="_")
pcst_isopodPS$Filter <- ifelse(grepl("*nohwe*",rownames(pcst_isopodPS)),
                             word(rownames(pcst_isopodPS),2,sep="_"),
                             pcst_isopodPS$Filter)
ggplot(pcst_isopodPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_isopodPS$Filter<-ifelse(pcst_isopodPS$Filter=="nohwe","No Filt",
                           ifelse(pcst_isopodPS$Filter=="out_any","HWE Out Any",
                                  ifelse(pcst_isopodPS$Filter=="out_all","HWE Out All","HWE Out Across")))
ggplot(pcst_isopodPS,aes(x=V1,y=Filter)) + theme_minimal() + ggtitle("specified bandwidth of 0.01") +
  xlab("PCst") +stat_density_ridges(bandwidth = 0.01)
ggplot(pcst_isopodPS[pcst_isopodPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("isopod_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_isopodPS$V1[pcst_isopodPS$Filter=="HWE Out All"],
       pcst_isopodPS$V1[pcst_isopodPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_isopodPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()


setwd("./vcfs/")
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
  t<-as.character(paste0("sim_",word(t[i],3,4,sep="_")))
  #a[i]<-(mean(as.numeric(as.matrix(fst_filesdat[[t]])),na.rm=T))
  b[i]<-mean(as.matrix(fst_simulations_dat[[i]][2:9]),na.rm=T)
  c[i]<-t
}
isopod_fsts<-data.frame(b,c)
colnames(isopod_fsts)<-c("Fst","Filter")
isopod_fsts$Filter<-ifelse(isopod_fsts$Filter=="sim_isopod_macfilt.recode","No Filt",
                             ifelse(isopod_fsts$Filter=="sim_out_any","HWE Out Any",
                                    ifelse(isopod_fsts$Filter=="sim_out_all","HWE Out All","HWE Out Across")))
ggplot(isopod_fsts,aes(x=Fst,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("isopod_PS_Fst") +
  xlab("Fst")
