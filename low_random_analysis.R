
setwd("/home/peawi142/HWE_Study_Dec2020/slim3.5_newsims/low_random_PS//")
low_random_PS_names <- list.files(path = "./randomized_pca_output/",pattern="*.csv",
                                   recursive = T,full.names = T)
low_random_PS_names<-low_random_PS_names[grepl("PCA",low_random_PS_names)]
low_random_PS <- lapply(low_random_PS_names,read.csv)
a<-paste(word(low_random_PS_names,4,sep="/"))
#a<-paste(word(a,1,2,3,sep="_"))
b<-paste(word(low_random_PS_names,6,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(low_random_PS) <- paste0(a,c)
low_random_PS<-low_random_PS[grepl("*subrep*",names(low_random_PS))]
dev.off();par(mfrow=c(1,5))
pop_list<-readLines("./randomized_pca_output/low_PS_6322479788840_/randomized/pop_map_randomized.csv")
plot(low_random_PS[[31]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="No Filt",xlim=c(-10,11),
     ylim=c(-6,6))
plot(low_random_PS[[11]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out All",xlim=c(-10,11),
     ylim=c(-6,6))
plot(low_random_PS[[21]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Any",xlim=c(-10,11),
     ylim=c(-6,6))
plot(low_random_PS[[1]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Across",xlim=c(-10,11),
     ylim=c(-6,6))

pcst_low_randomPS <- lapply(low_random_PS,pc_st,pop_list)
pcst_low_randomPS<-as.data.frame(do.call("rbind",pcst_low_randomPS))
pcst_low_randomPS$Filter<-rownames(pcst_low_randomPS)
pcst_low_randomPS$Filter <- word(pcst_low_randomPS$Filter,4,6,sep="_")
pcst_low_randomPS$Filter <- ifelse(grepl("*nohwe*",pcst_low_randomPS$Filter),
                                    word(pcst_low_randomPS$Filter,1,sep="_"),
                                    pcst_low_randomPS$Filter)
ggplot(pcst_low_randomPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_low_randomPS$Filter<-ifelse(pcst_low_randomPS$Filter=="nohwe","No Filt",
                                  ifelse(pcst_low_randomPS$Filter=="hwe_out_any","HWE Out Any",
                                         ifelse(pcst_low_randomPS$Filter=="hwe_out_all","HWE Out All","HWE Out Across")))
ggplot(pcst_low_randomPS,aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() + ggtitle("low_random_PS_PCA_PCst") +
  xlab("PCst")
ggplot(pcst_low_randomPS[pcst_low_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() + ggtitle("low_random_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_low_randomPS$V1[pcst_low_randomPS$Filter=="HWE Out All"],
       pcst_low_randomPS$V1[pcst_low_randomPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_low_randomPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()


fst_files<-list.files(path="../",pattern="*PS_fst*")
fst_files<-fst_files[grepl("low_PS",fst_files)]
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
low_randoms<-list.files(path="../rerun_christmas/")
low_randoms<-low_randoms[grepl("low_random",low_randoms)]
low_randoms<-low_randoms[grepl(".vcf",low_randoms)]
low_randoms<-paste0("../rerun_christmas/",low_randoms)
truefstsx<-lapply(truefsts_low,function(x){mean(as.matrix(x[2:7]),na.rm=T)})
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
ggplot(b,aes(x=Stat,y=Filt))+ geom_joy() + theme_minimal() + ggtitle("low_random_PS_PCA_PCst") +
  xlab("TrueFst-InfFst")
b$FileName<-names(fst_simulations_dat)
b$FileName<-word(b$FileName,1,2,sep="[.]")

pcst_low_randomPS$FileNames<-rownames(pcst_low_randomPS)
pcst_low_randomPS$FileNames<-word(pcst_low_randomPS$FileNames,1,2,sep="[.]")
pcst_low_randomPS$Fst<-b$InferredFst[match(pcst_low_randomPS$FileNames,b$FileName)]
pcst_low_randomPS$TrueFst<-b$TrueFst[match(pcst_low_randomPS$FileNames,b$FileName)]

dev.off()
ggplot(pcst_low_randomPS,aes(x=V1,y=Fst,col=Filter)) + geom_point()
library(cowplot)
low_random_fst_means <- ggplot(b,aes(x=InferredFst,y=Filt))+ geom_joy() + theme_minimal() +
  xlab("Inferred Fst")
low_random_fst_stand_truth <- ggplot(b,aes(x=Stat,y=Filt))+ geom_joy() + theme_minimal() +
  xlab("Inferred Fst")
low_random_pcst_fst <- ggplot(pcst_low_randomPS,aes(x=V1,y=Fst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
low_random_pcst_fst_noacross <- ggplot(pcst_low_randomPS[pcst_low_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=Fst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
low_random_pcst_all <- ggplot(pcst_low_randomPS,aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() +
  xlab("PCst")
low_random_pcst_noacross <- ggplot(pcst_low_randomPS[pcst_low_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_joy() + theme_minimal() +
  xlab("PCst")
low_random_pcst_truefst <- ggplot(pcst_low_randomPS,aes(x=V1,y=TrueFst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
low_random_pcst_truefst_noacross <- ggplot(pcst_low_randomPS[pcst_low_randomPS$Filter!="HWE Out Across",],aes(x=V1,y=TrueFst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
low_random_fst_truefst_noacross <- ggplot(pcst_low_randomPS[pcst_low_randomPS$Filter!="HWE Out Across",],aes(x=Fst,y=TrueFst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess')
low_random_fst_truefst <- ggplot(pcst_low_randomPS,aes(x=Fst,y=TrueFst,col=Filter)) + geom_point() + geom_smooth(method = 'loess')

plot_grid(low_random_fst_means, low_random_fst_stand_truth,
          low_random_pcst_fst,low_random_pcst_fst_noacross,
          low_random_pcst_all,low_random_pcst_noacross,low_random_pcst_truefst,
          low_random_fst_truefst, labels = c('InfFst_Means low_random PS', 'Fst_Standardized to Truth',
                                              'PCst vs InfFst', 'PCst vs InfFst Zoomed','PCst',
                                              'PCst zoomed','PCst vs TrueFst', 'Fst vs TrueFst'))
