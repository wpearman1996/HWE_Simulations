
setwd("/home/peawi142/HWE_Study_Dec2020/slim3.5_newsims/high_random_PS//")
high_random_PS_names <- list.files(path = "./randomized_pca_output/",pattern="*.csv",
                               recursive = T,full.names = T)
high_random_PS_names<-high_random_PS_names[grepl("PCA",high_random_PS_names)]
high_random_PS <- lapply(high_random_PS_names,read.csv)
a<-paste(word(high_random_PS_names,4,sep="/"))
#a<-paste(word(a,1,2,3,sep="_"))
b<-paste(word(high_random_PS_names,6,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(high_random_PS) <- paste0(a,c)
high_random_PS<-high_random_PS[grepl("*subrep*",names(high_random_PS))]
dev.off();par(mfrow=c(1,5))
pop_list<-readLines("./randomized_pca_output/high_PS_6322479788840_/randomized/pop_map_randomized.csv")
plot(high_random_PS[[31]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="No Filt",xlim=c(-10,11),
     ylim=c(-6,6))
plot(high_random_PS[[11]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out All",xlim=c(-10,11),
     ylim=c(-6,6))
plot(high_random_PS[[21]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Any",xlim=c(-10,11),
     ylim=c(-6,6))
plot(high_random_PS[[1]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Across",xlim=c(-10,11),
     ylim=c(-6,6))

pcst_high_randomPS <- lapply(high_random_PS,pc_st,pop_list)
pcst_high_randomPS<-as.data.frame(do.call("rbind",pcst_high_randomPS))
pcst_high_randomPS$Filter<-rownames(pcst_high_randomPS)
pcst_high_randomPS$Filter <- word(pcst_high_randomPS$Filter,4,6,sep="_")
pcst_high_randomPS$Filter <- ifelse(grepl("*nohwe*",pcst_high_randomPS$Filter),
                                word(pcst_high_randomPS$Filter,1,sep="_"),
                                pcst_high_randomPS$Filter)
ggplot(pcst_high_randomPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_high_randomPS$Filter<-ifelse(pcst_high_randomPS$Filter=="nohwe","No Filt",
                              ifelse(pcst_high_randomPS$Filter=="hwe_out_any","HWE Out Any",
                                     ifelse(pcst_high_randomPS$Filter=="hwe_out_all","HWE Out All","HWE Out Across")))
ggplot(pcst_high_randomPS,aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() + ggtitle("high_random_PS_PCA_PCst") +
  xlab("PCst")
ggplot(pcst_high_randomPS[pcst_high_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() + ggtitle("high_random_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_high_randomPS$V1[pcst_high_randomPS$Filter=="HWE Out All"],
       pcst_high_randomPS$V1[pcst_high_randomPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_high_randomPS, aes(x=V1, y=Filter)) + 
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
  file_gl_fst$pop<-as.factor(pop_list)
  fst_matrix<-stamppFst(file_gl_fst,nboots = 0)
  fst_matrix
}
high_randoms<-list.files(path="../rerun_christmas/")
high_randoms<-high_randoms[grepl("high_random",high_randoms)]
high_randoms<-high_randoms[grepl(".vcf",high_randoms)]
high_randoms<-paste0("../rerun_christmas/",high_randoms)
truefsts<-lapply(high_randoms,calc_fst_true)
names(truefsts)<-word(high_randoms,4,sep="_")
truefsts_high_random<-truefsts
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
ggplot(b,aes(x=Stat,y=Filt))+ geom_joy() + theme_minimal() + ggtitle("high_random_PS_PCA_PCst") +
  xlab("TrueFst-InfFst")
b$FileName<-names(fst_simulations_dat)
b$FileName<-word(b$FileName,1,2,sep="[.]")

pcst_high_randomPS$FileNames<-rownames(pcst_high_randomPS)
pcst_high_randomPS$FileNames<-word(pcst_high_randomPS$FileNames,1,2,sep="[.]")
pcst_high_randomPS$Fst<-b$InferredFst[match(pcst_high_randomPS$FileNames,b$FileName)]
pcst_high_randomPS$TrueFst<-b$TrueFst[match(pcst_high_randomPS$FileNames,b$FileName)]

dev.off()
ggplot(pcst_high_randomPS,aes(x=V1,y=Fst,col=Filter)) + geom_point()
library(cowplot)
high_random_fst_means <- ggplot(b,aes(x=InferredFst,y=Filt))+ geom_joy() + theme_minimal() +
  xlab("Inferred Fst")
high_random_fst_stand_truth <- ggplot(b,aes(x=Stat,y=Filt))+ geom_joy() + theme_minimal() +
  xlab("Inferred Fst")
high_random_pcst_fst <- ggplot(pcst_high_randomPS,aes(x=V1,y=Fst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
high_random_pcst_fst_noacross <- ggplot(pcst_high_randomPS[pcst_high_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=Fst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
high_random_pcst_all <- ggplot(pcst_high_randomPS,aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() +
  xlab("PCst")
high_random_pcst_noacross <- ggplot(pcst_high_randomPS[pcst_high_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() +
  xlab("PCst")
high_random_pcst_truefst <- ggplot(pcst_high_randomPS,aes(x=V1,y=TrueFst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
high_random_pcst_truefst_noacross <- ggplot(pcst_high_randomPS[pcst_high_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=TrueFst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
high_random_fst_truefst_noacross <- ggplot(pcst_high_randomPS[pcst_high_randomPS$Filter!="HWE Out Across",],aes(x=Fst,y=TrueFst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess')
high_random_fst_truefst <- ggplot(pcst_high_randomPS,aes(x=Fst,y=TrueFst,col=Filter)) + geom_point() + geom_smooth(method = 'loess')

plot_grid(high_random_fst_means, high_random_fst_stand_truth,
          high_random_pcst_fst,high_random_pcst_fst_noacross,
          high_random_pcst_all,high_random_pcst_noacross,high_random_pcst_truefst,
          high_random_fst_truefst, labels = c('InfFst_Means high_random PS', 'Fst_Standardized to Truth',
                                          'PCst vs InfFst', 'PCst vs InfFst Zoomed','PCst',
                                          'PCst zoomed','PCst vs TrueFst', 'Fst vs TrueFst'))
