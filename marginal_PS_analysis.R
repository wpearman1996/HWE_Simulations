
setwd("/home/peawi142/HWE_Simulations/heterozgosity_stats/marg_PS/")
marg_PS_names <- list.files(path = "./het_calcs_pca_output/",pattern="*.csv",
                           recursive = T,full.names = T)
marg_PS_names<-marg_PS_names[grepl("PCA",marg_PS_names)]
marg_PS <- lapply(marg_PS_names,read.csv)
a<-paste(word(marg_PS_names,4,sep="/"))
#a<-paste(word(a,1,2,3,sep="_"))
b<-paste(word(marg_PS_names,6,sep="/"))
c<-paste(word(b,1,2,sep="[.]"))
names(marg_PS) <- paste0(a,c)
marg_PS<-marg_PS[grepl("*subrep*",names(marg_PS))]
dev.off();par(mfrow=c(1,5))
pop_list<-c(rep("p1",30),rep("p2",30),rep("p3",30),rep("p4",30),rep("p5",30),rep("p6",30))
plot(marg_PS[[31]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="No Filt",xlim=c(-10,11),
     ylim=c(-6,6))
plot(marg_PS[[11]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out All",xlim=c(-10,11),
     ylim=c(-6,6))
plot(marg_PS[[21]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Any",xlim=c(-10,11),
     ylim=c(-6,6))
plot(marg_PS[[1]][2:3],col=as.factor(pop_list),pch=19,cex=2,main="HWE Out Across",xlim=c(-10,11),
     ylim=c(-6,6))

pcst_margPS <- lapply(marg_PS,pc_st,pop_list)
pcst_margPS<-as.data.frame(do.call("rbind",pcst_margPS))
pcst_margPS$Filter<-rownames(pcst_margPS)
pcst_margPS$Filter <- word(pcst_margPS$Filter,4,6,sep="_")
pcst_margPS$Filter <- ifelse(grepl("*nohwe*",pcst_margPS$Filter),
                            word(pcst_margPS$Filter,1,sep="_"),
                            pcst_margPS$Filter)
ggplot(pcst_margPS, aes(x=V1, y=Filter)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0))  + coord_flip()

pcst_margPS$Filter<-ifelse(pcst_margPS$Filter=="nohwe","No Filt",
                          ifelse(pcst_margPS$Filter=="hwe_out_any","HWE Out Any",
                                 ifelse(pcst_margPS$Filter=="hwe_out_all","HWE Out All","HWE Out Across")))
ggplot(pcst_margPS,aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("marg_PS_PCA_PCst") +
  xlab("PCst")
ggplot(pcst_margPS,aes(x=V1,y=Filter)) + theme_minimal() + ggtitle("marg_PS_PCA_PCst") +
  xlab("PCst") + stat_density_ridges(bandwidth=0.003)
ggplot(pcst_margPS[pcst_margPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() + ggtitle("marg_PS_PCA_PCst") +
  xlab("PCst")
t.test(pcst_margPS$V1[pcst_margPS$Filter=="HWE Out All"],
       pcst_margPS$V1[pcst_margPS$Filter=="No Filt"],alternative = "greater")
ggplot(pcst_margPS, aes(x=V1, y=Filter)) + 
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
#  c[i]<-t
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
marg_fsts<-b
pcst_margPS$FileNames<-rownames(pcst_margPS)
pcst_margPS$FileNames<-word(pcst_margPS$FileNames,1,2,sep="[.]")
pcst_margPS$Fst<-b$InferredFst[match(pcst_margPS$FileNames,b$FileName)]

dev.off()
ggplot(pcst_margPS,aes(x=V1,y=Fst,col=Filter)) + geom_point()
library(cowplot)
marg_fst_means <- ggplot(b,aes(x=InferredFst,y=Filt))+ geom_density_ridges() + theme_minimal() +
  xlab("Inferred Fst")
marg_pcst_fst <- ggplot(pcst_margPS,aes(x=V1,y=Fst,col=Filter)) + geom_point() + geom_smooth(method = 'loess') +
  xlab("PCst")
marg_pcst_fst_noacross <- ggplot(pcst_margPS[pcst_margPS$Filter!="HWE Out Across",],aes(x=V1,y=Fst,col=Filter)) +
  geom_point() + geom_smooth(method = 'loess') + xlab("PCst")
marg_pcst_all <- ggplot(pcst_margPS,aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() +
  xlab("PCst")
marg_pcst_noacross <- ggplot(pcst_margPS[pcst_margPS$Filter!="HWE Out Across",],aes(x=V1,y=Filter))+ geom_density_ridges() + theme_minimal() +
  xlab("PCst")
pcst_margPS$Rep<-word(pcst_margPS$FileNames,1,3,sep="_")
p <- ggplot(pcst_margPS, aes(x=Filter, y=Fst, group=Rep)) +
  geom_line(aes(color=Rep))+
  geom_point(aes(color=Rep)) + scale_color_brewer(palette="Paired")+
  theme_minimal()
p
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
df3 <- data_summary(pcst_margPS, varname="Fst", 
                    groupnames=c("Filter", "Rep"))

ggplot(df3, aes(x=Filter, y=Fst, group=Rep, color=Rep)) + 
  geom_errorbar(aes(ymin=Fst-sd, ymax=Fst+sd), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()


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
out_across<-admix_plot(clummped_marg_dat$out_across,10,180,6,F,cbbPalette[c(1:5,7)],"Out Across")
out_any<-admix_plot(clummped_marg_dat$out_any,10,180,6,F,cbbPalette[c(1:5,7)],"Out Any")
out_all<-admix_plot(clummped_marg_dat$out_all,10,180,6,F,cbbPalette[c(1:5,7)],"Out All")
nofilt<-admix_plot(clummped_marg_dat$out_any,10,180,6,F,cbbPalette[c(1:5,7)],"No Filt")
clummped_marg_dat_melt<-do.call("rbind",clummped_marg_dat)
clummped_marg_dat_melt$Filter<-word(rownames(clummped_marg_dat_melt),1,sep="[.]")
meltq <- reshape::melt(clummped_marg_dat_melt,id.vars=c("Filter"))
meltq$IndName<-rep(1:180,4)
meltq$Filter = factor(meltq$Filter, levels=c('nofilt','out_any','out_all','out_across'))
marg_struc_nucdist <- ggplot(struc_mats,aes(x=struc_mats,y=filter)) + theme_minimal() + xlab("Nuc Dist") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1)
marg_struc_mats <- struc_mats
marg_struc_nucdist
library(cowplot)
plot_grid(marg_struc_nucdist,marg_struc_nucdist)
clummped_marg_dat <- run_structure_analysis("./", k=6, pop_list=pop_list,
                                           simulation="marg_PS",useclumpp=F)
Struc_plot_marg<-ggplot(data=meltq,aes(x=IndName,y=value,fill=variable))+
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
#Struc_plot_marg
#marg_struc_nucdist
ggdraw()+
  draw_plot(marg_struc_nucdist,0,0,.4,1)+
  draw_plot(Struc_plot_marg + theme(panel.spacing = unit(1.5, "lines")),0.55,.055,0.55,0.9,hjust = .3,vjust = -.11)

plot_grid(marg_fst_means, marg_pcst_fst,marg_pcst_fst_noacross,
          marg_pcst_all,marg_pcst_noacross,marg_struc_nucdist
          , labels = c('InfFst_Means marg PS','PCst vs InfFst', 'PCst vs InfFst Zoomed','PCst',
                       'PCst zoomed', 'Structure Nuc Dist'))


#### Heterozygosit and basic stats
setwd("/home/peawi142/HWE_Simulations/heterozgosity_stats/marg_PS/")
files_stat_names <- list.files(path = "./het_calcs_pca_output/",pattern="*summary_stats.csv",
                               recursive = T,full.names = T)
files_stat_names<-files_stat_names[grepl("subrep",files_stat_names)]
files_stat<-lapply(files_stat_names,read.csv)
extract_fis<-function(file){
  file$x
}
marg_fis_dat<-lapply(files_stat,extract_fis)
newnames<-word(files_stat_names,6,sep="/") %>%
  word(.,2,3,sep="_") 
marg_fis_dat<-data.frame(do.call("rbind",marg_fis_dat));colnames(marg_fis_dat)<-c("Ho","Hs","Ht","Dst","Dstp","Fst","Fstp","Fis","Dest")
marg_fis_dat$newname<-newnames
marg_fis_dat$Filt<-ifelse(marg_fis_dat$newname=="out_across","HWE Out Across",ifelse(marg_fis_dat$newname=="out_any","HWE Out Any",
                                                                                     ifelse(marg_fis_dat$newname=="out_all","HWE Out All", "No HWE")))

marg_fis_distr <- ggplot(marg_fis_dat,aes(x=Ho,y=Filt)) + theme_minimal() + xlab("Fis") + ylab("HWE Filter") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), scale=1)
marg_fis_distr
