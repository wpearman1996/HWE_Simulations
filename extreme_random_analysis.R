
setwd("/home/peawi142/HWE_Simulations/randomized/extreme_rand//")
extreme_rand_PS_names <- list.files(path = "./randomized_pca_output/",pattern="*.csv",
                                recursive = T,full.names = T)
extreme_rand_PS_names<-extreme_rand_PS_names[grepl("PCA",extreme_rand_PS_names)]
extreme_rand_PS <- lapply(extreme_rand_PS_names,read.csv)
a<-paste(word(extreme_rand_PS_names,4,sep="/"))
#a<-paste(word(a,1,2,3,sep="_"))
b<-paste(word(extreme_rand_PS_names,6,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(extreme_rand_PS) <- paste0(a,c)
extreme_rand_PS<-extreme_rand_PS[grepl("*subrep*",names(extreme_rand_PS))]
dev.off();par(mfrow=c(1,5))
pop_list<-readLines("./randomized_pca_output/extreme_PS_6322559800512_/randomized/pop_map_randomized.csv")
plot(extreme_rand_PS[[31]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="No Filt",xlim=c(-10,11),
     ylim=c(-6,6))
plot(extreme_rand_PS[[11]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out All",xlim=c(-10,11),
     ylim=c(-6,6))
plot(extreme_rand_PS[[21]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Any",xlim=c(-10,11),
     ylim=c(-6,6))
plot(extreme_rand_PS[[1]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Across",xlim=c(-10,11),
     ylim=c(-6,6))

pcst_extreme_randPS <- lapply(extreme_rand_PS,pc_st,pop_list)
pcst_extreme_randPS<-as.data.frame(do.call("rbind",pcst_extreme_randPS))
pcst_extreme_randPS$Filter<-rownames(pcst_extreme_randPS)
pcst_extreme_randPS$Filter <- word(pcst_extreme_randPS$Filter,4,6,sep="_")
pcst_extreme_randPS$Filter <- ifelse(grepl("*nohwe*",pcst_extreme_randPS$Filter),
                                 word(pcst_extreme_randPS$Filter,1,sep="_"),
                                 pcst_extreme_randPS$Filter)
ggplot(pcst_extreme_randPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_extreme_randPS$Filter<-ifelse(pcst_extreme_randPS$Filter=="nohwe","No Filt",
                               ifelse(pcst_extreme_randPS$Filter=="hwe_out_any","HWE Out Any",
                                      ifelse(pcst_extreme_randPS$Filter=="hwe_out_all","HWE Out All","HWE Out Across")))
ggplot(pcst_extreme_randPS,aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("extreme_rand_PS_PCA_PCst") +
  xlab("PCst")
ggplot(pcst_extreme_randPS,aes(x=V1,y=Filter)) + theme_minimal() + ggtitle("extreme_rand_PS_PCA_PCst") +
  xlab("PCst") + stat_density_ridges(bandwidth=0.003)
ggplot(pcst_extreme_randPS[pcst_extreme_randPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("extreme_rand_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_extreme_randPS$V1[pcst_extreme_randPS$Filter=="HWE Out All"],
       pcst_extreme_randPS$V1[pcst_extreme_randPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_extreme_randPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

fst_simulations<-list.files(pattern="*csv",recursive = T,full.names = T)
fst_simulations <- fst_simulations[grepl("Fst",fst_simulations)]
fst_simulations_dat<-lapply(fst_simulations,read.csv)
names(fst_simulations_dat) <- paste0(word(fst_simulations,3,sep="/"),word(fst_simulations,5,sep="/"))
calc_dist_mean <- function(true_fst,inferred_fst){
  mean(true_fst)-mean(inferred_fst)
}
fst_simulations_dat<-fst_simulations_dat[grepl("subrep",names(fst_simulations_dat))]
xtest<-list()
b<-vector()
for ( i in 1:length(fst_simulations_dat)){
  b[i]<-mean(as.matrix(fst_simulations_dat[[i]][2:7]),na.rm=T)
  c[i]<-t
}

b<-data.frame(b)
colnames(b)<-"InferredFst"
b$Name<-names(fst_simulations_dat)
b$Seed<-word(b$Name,3,sep="_")
b$Filt<-word(b$Name,5,6,sep="_")
b$Filt<-ifelse(b$Filt=="out_across","HWE Out Across",ifelse(b$Filt=="out_any","HWE Out Any",
                                                            ifelse(b$Filt=="out_all","HWE Out All", "No HWE")))
b$FileName<-names(fst_simulations_dat)
b$FileName<-word(b$FileName,1,2,sep="[.]")
extreme_rand_Fst<-b
pcst_extreme_randPS$FileNames<-rownames(pcst_extreme_randPS)
pcst_extreme_randPS$FileNames<-word(pcst_extreme_randPS$FileNames,1,2,sep="[.]")
pcst_extreme_randPS$Fst<-b$InferredFst[match(pcst_extreme_randPS$FileNames,b$FileName)]

dev.off()
ggplot(pcst_extreme_randPS,aes(x=V1,y=Fst,col=Filter)) + geom_point()
library(cowplot)
extreme_rand_fst_means <- ggplot(b,aes(x=InferredFst,y=Filt))+ geom_density_ridges() + theme_minimal() +
  xlab("Inferred Fst")
extreme_rand_pcst_fst <- ggplot(pcst_extreme_randPS,aes(x=V1,y=Fst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
extreme_rand_pcst_fst_noacross <- ggplot(pcst_extreme_randPS[pcst_extreme_randPS$Filter!="HWE Out Across",],aes(x=V1,y=Fst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
extreme_rand_pcst_all <- ggplot(pcst_extreme_randPS,aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() +
  xlab("PCst")
extreme_rand_pcst_noacross <- ggplot(pcst_extreme_randPS[pcst_extreme_randPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() +
  xlab("PCst")
####### STRUCTURE
setwd("./structure_output/")
files<-list.files(pattern="*out_f")
struc_mats<-lapply(files,read_struc_nucmat)
struc_mats<-do.call("rbind",struc_mats)
struc_mats<-data.frame(struc_mats,files)
struc_mats$filter<-ifelse(grepl("nohwe",struc_mats$files),"No Filter",
                          ifelse(grepl("out_any",struc_mats$files),"HWE Out Any",
                                 ifelse(grepl("out_across",struc_mats$files),
                                        "HWE Out Across", "HWE Out All"))) 

dev.off();par(mfrow=c(1,4))
out_across<-admix_plot(clummped_extreme_rand_dat$out_across,10,180,6,F,cbbPalette[c(1:5,7)],"Out Across")
out_any<-admix_plot(clummped_extreme_rand_dat$out_any,10,180,6,F,cbbPalette[c(1:5,7)],"Out Any")
out_all<-admix_plot(clummped_extreme_rand_dat$out_all,10,180,6,F,cbbPalette[c(1:5,7)],"Out All")
nofilt<-admix_plot(clummped_extreme_rand_dat$out_any,10,180,6,F,cbbPalette[c(1:5,7)],"No Filt")
clummped_extreme_rand_dat_melt<-do.call("rbind",clummped_extreme_rand_dat)
clummped_extreme_rand_dat_melt$Filter<-word(rownames(clummped_extreme_rand_dat_melt),1,sep="[.]")
meltq <- reshape::melt(clummped_extreme_rand_dat_melt,id.vars=c("Filter"))
meltq$IndName<-rep(1:180,4)
meltq$Filter = factor(meltq$Filter, levels=c('nofilt','out_any','out_all','out_across'))
extreme_rand_struc_nucdist <- ggplot(struc_mats,aes(x=struc_mats,y=filter)) + theme_minimal() + xlab("Nuc Dist") + ylab("HWE Filter") #+
extreme_rand_struc_mats <- struc_mats

extreme_rand_struc_nucdist
library(cowplot)
plot_grid(extreme_rand_struc_nucdist,extreme_rand_struc_nucdist)
clummped_extreme_rand_dat <- run_structure_analysis("./", k=6, pop_list=pop_list,
                                                simulation="extreme_rand_PS",useclumpp=F)
Struc_plot_extreme_rand<-ggplot(data=meltq,aes(x=IndName,y=value,fill=variable))+
  facet_wrap(Filter~.,nrow=4,strip.position = NULL)+
  ggtitle("  ")+
  theme_classic()+theme(axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        strip.background = element_blank(),
                        axis.text.x=element_blank(),
                        axis.text = element_blank(),
                        strip.text.y.left = element_text(angle = 0),
                        rect = element_blank(),axis.line.x = element_blank(),
                        axis.line.y = element_blank(),
                        legend.position = "none", strip.text.x = element_blank() )+
  ylab("")+xlab("")+
  scale_fill_manual(values = cbp1[c(1:6)])+
  geom_bar(stat="identity",width=1,col="black",lwd=0.00)
#Struc_plot_extreme_rand
#extreme_rand_struc_nucdist
ggdraw()+
  draw_plot(extreme_rand_struc_nucdist,0,0,.4,1)+
  draw_plot(Struc_plot_extreme_rand + theme(panel.spacing = unit(1.5, "lines")),0.55,.055,0.55,0.9,hjust = .3,vjust = -.11)

plot_grid(extreme_rand_fst_means, extreme_rand_pcst_fst,extreme_rand_pcst_fst_noacross,
          extreme_rand_pcst_all,extreme_rand_pcst_noacross,extreme_rand_struc_nucdist
          , labels = c('InfFst_Means extreme_rand PS','PCst vs InfFst', 'PCst vs InfFst Zoomed','PCst',
                       'PCst zoomed', 'Structure Nuc Dist'))


#### Heterozygosit and basic stats


