
setwd("/home/peawi142/HWE_Study_Dec2020/slim3.5_newsims/high_PS/")
high_PS_names <- list.files(path = "./pca_output/",pattern="*.csv",
                               recursive = T,full.names = T)
high_PS_names<-high_PS_names[grepl("PCA",high_PS_names)]
high_PS <- lapply(high_PS_names,read.csv)
a<-paste(word(high_PS_names,4,sep="/"))
#a<-paste(word(a,1,2,3,sep="_"))
b<-paste(word(high_PS_names,6,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(high_PS) <- paste0(a,c)
high_PS<-high_PS[grepl("*subrep*",names(high_PS))]
dev.off();par(mfrow=c(1,5))
pop_list<-c(rep("p1",30),rep("p2",30),rep("p3",30),rep("p4",30),rep("p5",30),rep("p6",30))
plot(high_PS[[31]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="No Filt",xlim=c(-10,11),
     ylim=c(-6,6))
plot(high_PS[[11]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out All",xlim=c(-10,11),
     ylim=c(-6,6))
plot(high_PS[[21]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Any",xlim=c(-10,11),
     ylim=c(-6,6))
plot(high_PS[[1]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Across",xlim=c(-10,11),
     ylim=c(-6,6))

pcst_highPS <- lapply(high_PS,pc_st,pop_list)
pcst_highPS<-as.data.frame(do.call("rbind",pcst_highPS))
pcst_highPS$Filter<-rownames(pcst_highPS)
pcst_highPS$Filter <- word(pcst_highPS$Filter,4,6,sep="_")
pcst_highPS$Filter <- ifelse(grepl("*nohwe*",pcst_highPS$Filter),
                                word(pcst_highPS$Filter,1,sep="_"),
                                pcst_highPS$Filter)
ggplot(pcst_highPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_highPS$Filter<-ifelse(pcst_highPS$Filter=="nohwe","No Filt",
                              ifelse(pcst_highPS$Filter=="hwe_out_any","HWE Out Any",
                                     ifelse(pcst_highPS$Filter=="hwe_out_all","HWE Out All","HWE Out Across")))
ggplot(pcst_highPS,aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() + ggtitle("high_PS_PCA_PCst") +
  xlab("PCst")
ggplot(pcst_highPS[pcst_highPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() + ggtitle("high_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_highPS$V1[pcst_highPS$Filter=="HWE Out All"],
       pcst_highPS$V1[pcst_highPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_highPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()


fst_files<-list.files(path="../",pattern="*PS_fst*")
fst_files<-fst_files[grepl("high_PS",fst_files)]
fst_filesdat<-lapply(paste0("../",fst_files),readLines)
names(fst_filesdat)<-paste0("sim_",word(fst_files,5,sep="_"))
make_fst_matrix <- function(vector){
  x<-matrix(NA,nrow=6,ncol=6)
  colnames(x)<-c("P1","P3","P5","P7","P9","P11")
  rownames(x)<-c("P1","P3","P5","P7","P9","P11")
  vector<-vector[16:30]
  vector<-as.numeric(vector)
  x[2,1]<-vector[1];x[3,1]<-vector[2];x[4,1]<-vector[3];x[5,1]<-vector[4];x[6,1]<-vector[5]
  x[3,2]<-vector[6];x[4,2]<-vector[7];x[5,2]<-vector[8];x[6,2]<-vector[9]
  x[4,3]<-vector[10];x[5,3]<-vector[11];x[6,3]<-vector[12]
  x[5,4]<-vector[13];x[6,4]<-vector[14]
  x[6,5]<-vector[15]
  x
}
fst_filesdat<-lapply(fst_filesdat,make_fst_matrix)
fst_simulations<-list.files(pattern="*csv",recursive = T,full.names = T)
fst_simulations <- fst_simulations[grepl("Fst",fst_simulations)]
fst_simulations_dat<-lapply(fst_simulations,read.csv)
names(fst_simulations_dat) <- paste0(word(fst_simulations,3,sep="/"),word(fst_simulations,5,sep="/"))
calc_dist_mean <- function(true_fst,inferred_fst){
  mean(true_fst)-mean(inferred_fst)
}
fst_simulations_dat<-fst_simulations_dat[grepl("subrep",names(fst_simulations_dat))]
xtest<-list()
a<-vector()
b<-vector()
c<-vector()
for ( i in 1:length(fst_simulations_dat)){
  t<-names(fst_simulations_dat)
  t<-as.character(paste0("sim_",word(t[i],3,sep="_")))
  a[i]<-(mean(as.numeric(as.matrix(fst_filesdat[[t]])),na.rm=T))
  b[i]<-mean(as.matrix(fst_simulations_dat[[i]][2:7]),na.rm=T)
  c[i]<-t
}
names(a)<-c
a<-data.frame(a)
colnames(a)<-"TrueFst"
a$Seed<-word(c,2,sep="_")
calc_fst_true<-function(vcf){
  file <- read.vcfR(vcf,verbose = FALSE)
  file_gl <- vcfR2genlight(file)
  file_gl_fst<-file_gl
  file_gl_fst$pop<-as.factor(c(rep("p1",30),rep("p2",30),
                               rep("p3",30),rep("p4",30),
                               rep("p5",30),rep("p6",30)))
  fst_matrix<-stamppFst(file_gl_fst,nboots = 0)
  fst_matrix
}
highs<-list.files(path="../rerun_christmas/")
highs<-highs[grepl("high",highs)]
highs<-highs[grepl(".vcf",highs)]
highs<-paste0("../rerun_christmas/",highs)
truefsts<-lapply(highs,calc_fst_true)
names(truefsts)<-word(highs,4,sep="_")
truefsts_high<-truefsts
truefstsx<-lapply(truefsts_high,function(x){mean(as.matrix(x[2:7]),na.rm=T)})
truefstsx<-do.call("rbind",truefstsx)
truefstsx<-data.frame(truefstsx)
truefstsx$seed<-rownames(truefstsx)
a<-unique(a)
b<-data.frame(b)
colnames(b)<-"InferredFst"
b$Name<-names(fst_simulations_dat)
b$Seed<-word(b$Name,3,sep="_")
b$Filt<-word(b$Name,5,6,sep="_")
b$Filt<-ifelse(b$Filt=="out_across","HWE Out Across",ifelse(b$Filt=="out_any","HWE Out Any",
                                                            ifelse(b$Filt=="out_all","HWE Out All", "No HWE")))
b$TrueFst<-truefstsx$truefstsx[match(b$Seed,truefstsx$seed)]
b$Stat <- b$TrueFst-b$InferredFst
ggplot(b,aes(x=Stat,y=Filt))+ geom_joy() + theme_minimal() + ggtitle("high_PS_PCA_PCst") +
  xlab("TrueFst-InfFst")
b$FileName<-names(fst_simulations_dat)
b$FileName<-word(b$FileName,1,2,sep="[.]")

pcst_highPS$FileNames<-rownames(pcst_highPS)
pcst_highPS$FileNames<-word(pcst_highPS$FileNames,1,2,sep="[.]")
pcst_highPS$Fst<-b$InferredFst[match(pcst_highPS$FileNames,b$FileName)]
pcst_highPS$TrueFst<-b$TrueFst[match(pcst_highPS$FileNames,b$FileName)]

dev.off()
ggplot(pcst_highPS,aes(x=V1,y=Fst,col=Filter)) + geom_point()
library(cowplot)
high_fst_means <- ggplot(b,aes(x=InferredFst,y=Filt))+ geom_joy() + theme_minimal() +
  xlab("Inferred Fst")
high_fst_stand_truth <- ggplot(b,aes(x=Stat,y=Filt))+ geom_joy() + theme_minimal() +
  xlab("Inferred Fst")
high_pcst_fst <- ggplot(pcst_highPS,aes(x=V1,y=Fst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
high_pcst_fst_noacross <- ggplot(pcst_highPS[pcst_highPS$Filter!="HWE Out Across",],aes(x=V1,y=Fst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
high_pcst_all <- ggplot(pcst_highPS,aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() +
  xlab("PCst")
high_pcst_noacross <- ggplot(pcst_highPS[pcst_highPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() +
  xlab("PCst")
high_pcst_truefst <- ggplot(pcst_highPS,aes(x=V1,y=TrueFst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
high_pcst_truefst_noacross <- ggplot(pcst_highPS[pcst_highPS$Filter!="HWE Out Across",],aes(x=V1,y=TrueFst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
high_fst_truefst_noacross <- ggplot(pcst_highPS[pcst_highPS$Filter!="HWE Out Across",],aes(x=Fst,y=TrueFst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess')
high_fst_truefst <- ggplot(pcst_highPS,aes(x=Fst,y=TrueFst,col=Filter)) + geom_point() + geom_smooth(method = 'loess')

plot_grid(high_fst_means, high_fst_stand_truth,
          high_pcst_fst,high_pcst_fst_noacross,
          high_pcst_all,high_pcst_noacross,high_pcst_truefst,
          high_fst_truefst, labels = c('InfFst_Means high PS', 'Fst_Standardized to Truth',
                                                       'PCst vs InfFst', 'PCst vs InfFst Zoomed','PCst',
                                                       'PCst zoomed','PCst vs TrueFst', 'Fst vs TrueFst'))