read_struc_nucmat<-function(file){
  x<-readLines(file)
  x<-x[36:42]
  data<-list()
  for (i in 3:length(x) ) {
    data[[i]]<-read.table(textConnection(x[[i]]))
  }
  mat<-as.numeric(as.matrix(do.call("rbind",data)[,2:7]))
  mean(na.omit(mat))
}

library(hierfstat)
library(tidyverse)
library(adegenet)
library(dartR)
library(pegas)
library(vcfR)
library(dartR)
library(RColorBrewer)
library(pophelper)
admix_plot<-function(clump_align_tab,nrep,nind,title){
  line_locs<-c(30,60,90,120,150)
  clump_align_tab$V1<-NULL
  clump_align_tab[,ncol(clump_align_tab)]<-NULL
  #x<-split(clump_align_tab,sort(rep(seq(1:nrep),nind)))
  #x<-lapply(x,as.matrix)
  #z<-Reduce("+", x) / length(x)
  z<-clump_align_tab
  #col_pal<-col+pal
  col_pal<-brewer.pal(6,"Paired")
  #cols<-col_pal#[1:ncol(z)]
  barplot(t(z),col=col_pal,xlab="Admixture Proportions",ylab=NULL,
          space= 0,border=NA,axisnames = F,horiz=T, main = title,
          cex.axis = 1.2, cex.names = 1.2,cex.lab=1.2,cex.main=1.5);abline(h=line_locs)
}

clumppExport <- function(qlist=NULL,prefix=NA,parammode=NA,paramrep=NA,useexe=FALSE,dir_name)
{
  # check input
  is.qlist(qlist)
  
  if(is.na(prefix)) prefix <- "pop"
  prefix <- paste0(prefix,"_K")
  if(!is.logical(useexe)) stop("clumppExport: Argument 'useexe' set incorrectly. Set as TRUE or FALSE.")
  
  # get tabulated runs
  df1 <- pophelper::tabulateQ(qlist)
  df2 <- pophelper::summariseQ(df1)
  df1l <- as.list(df1)
  df2l <- as.list(df2)
  
  if(is.null(names(qlist))) names(qlist) <- paste0("sample",1:length(qlist))
  
  # k val duplicated
  if(any(duplicated(df2l$k))) stop("clumppExport: Repeating values of K found.")
  # do ind vary?
  if(!all(df2l$ind[1]==df2l$ind)) warning("clumppExport: Number of individuals vary between runs.")
  
  e <- 1
  p <- 1
  len1 <- length(df2l$k)
  while (e <= len1)
  {
    k <- df2l$k[e]
    ind <- df2l$ind[e]
    runs <- df2l$runs[e]
    
    ldata <- vector("list",length=runs)
    for (f in 1:runs)
    {
      sel <- which(names(qlist)==as.character(df1l$file[p]))
      dframe1 <- qlist[[sel]]
      
      # generate df
      dframe3 <- as.matrix(data.frame(V1=paste0(1:ind,":"),dframe1,last=as.character(rep(1,ind)),stringsAsFactors=FALSE))
      
      # add dataframes to list
      ldata[[f]] <- dframe3
      rm(dframe3)
      p=p+1
    }
    
    if(runs > 1 && k > 1)
    {
      currwd <- getwd()
      if(as.numeric(file.access(currwd,2))==-1) stop(paste0("clumppExport: Directory ",currwd," has no write permission."))
      
      dir.create(paste0(currwd,"/",dir_name))
      setwd(paste0(currwd,"/",dir_name))
      cat(paste0("Folder created: ",basename(getwd()),"\n"))  
      out <- paste0(prefix,k,"-combined.txt")
      
      ## file output block
      
      # make 2 line space
      spacer <- matrix(rep("  ",(k+2)*2),nrow=2)
      
      # write file
      write(t(format(ldata[[1]],nsmall=15)),paste(out),ncolumns=k+2)
      for (i in 2:length(ldata))
      {
        write(t(spacer),paste(out),ncolumns=k+2,append=TRUE)
        write(t(format(ldata[[i]],nsmall=15)),append=TRUE,paste(out),ncolumns=k+2)
      }
      cat(paste0(out),"exported.\n")
      
      ## paramfile section
      T1 <- factorial(k)*((length(ldata)*(length(ldata)-1))/2)*k*ind
      if(T1 <= 100000000)
      {
        if(is.na(parammode)) parammode <- 2
        if(is.na(paramrep)) paramrep <- 20
      }else{
        if(is.na(parammode)) parammode <- 3
        if(is.na(paramrep)) paramrep <- 500
      }
      out1 <- base::gsub(".txt","",out)
      params <- c("DATATYPE 1 ",
                  "INDFILE NOTNEEDED.indfile ",
                  paste0("POPFILE ",out," "),
                  paste0("OUTFILE ",out1,"-merged.txt "),
                  paste0("MISCFILE ",out1,"-miscfile.txt "),
                  paste0("K ",k," "),
                  paste0("C ",ind," "),
                  paste0("R ",length(ldata)," "),
                  paste0("M ",parammode," "),
                  "W 0 ",
                  "S 2 ",
                  "GREEDY_OPTION 2 ",
                  paste0("REPEATS ",paramrep," "),
                  "PERMUTATIONFILE NOTNEEDED.permutationfile ",
                  "PRINT_PERMUTED_DATA 1 ",
                  paste0("PERMUTED_DATAFILE ",out1,"-aligned.txt "),
                  "PRINT_EVERY_PERM 0 ",
                  paste0("EVERY_PERMFILE ",out1,".every_permfile "),
                  "PRINT_RANDOM_INPUTORDER 0 ",
                  paste0("RANDOM_INPUTORDERFILE ",out1,".random_inputorderfile "),
                  "OVERRIDE_WARNINGS 0 ",
                  "ORDER_BY_RUN 0 ")
      
      write(params,"paramfile")
      cat(paste0("paramfile exported.\n"))
      
      # autorun clumpp executable
      if(useexe)
      {
        # identify OS
        sysos <- "unix64"
        if(sysos=="unix64")
        {
          #file.copy(system.file("/bin/CLUMPP",package="pophelper"),".")
          #system("chmod 777 CLUMPP")
          system("CLUMPP")
          #unlink("CLUMPP",force=TRUE)
        }
        
        # if OS is unidentified, give error
        if(sysos=="unknown") warning("clumppExport: CLUMPP executable not run because system cannot be identified as windows, mac or linux.")
      }
      
      setwd(paste(currwd))
      cat("-----------------------\n")
    }else
    {
      if(k==1) message(paste0(prefix,k," not exported. K less than 2.\n"))
      if(runs < 2) message(paste0(prefix,k," not exported. Repeats less than 2.\n"))
      cat("-----------------------\n")
    }
    e <- e + 1
  }
  
  cat("Run completed.\n")
}

genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

gen_rand <- function(genind_tab,nloci){
  x<-1:ncol(genind_tab)
  x<-x[c(TRUE,FALSE)]
  x<-sample(x,nloci)
  x<-c(x,x+1)
  x<-x[order(x)]
  y<-genind_tab[,x]
  y
}
check_hw_per_pop <- function(genlight,population){
  y <- as.character(unique(genlight$pop)) 
  y <- y[! y %in% population]
  x<-gl.drop.pop(genlight,y,v=0,mono.rm = T)
  t <- gl2gi(x)
  data_hw<-hw.test(t,B=0)
  data_hw<-data.frame(data_hw)
  data_hw_sig<-data_hw[p.adjust(data_hw$Pr.chi.2...,"BH") <= 0.05,]
  #  data_hw_sig<-data_hw[data_hw$Pr.chi.2... <= 0.05/x$n.loc,]
  rownames(data_hw_sig)
}
check_hw_per_pop_across <- function(genlight){
  # y <- as.character(unique(genlight$pop)) 
  #  y <- y[! y %in% population]
  #  x<-gl.drop.pop(genlight,y,v=0)
  t <- gl2gi(genlight)
  t$pop <- as.factor(rep("p1",nrow(t$tab)))
  data_hw<-hw.test(t,B=0)
  data_hw<-data.frame(data_hw)
  data_hw_sig<-data_hw[p.adjust(data_hw$Pr.chi.2...,"BH") <= 0.05,]
  #  data_hw_sig<-data_hw[data_hw$Pr.chi.2... <= 0.05/x$n.loc,]
  rownames(data_hw_sig)
}
pc_st <- function(pc_matrix,pop_list){
  pc_matrix[,1:5]<-apply(X = pc_matrix[,1:5],FUN = as.numeric,MARGIN=2)
  pc_matrix<-pc_matrix[,1:5]
  pc_matrix<-data.frame(pc_matrix)
  pc_matrix[,6]<-pop_list
  out <- c()
  for(k in unique(pc_matrix[,6])){
    within_dist <- mean(dist(pc_matrix[pc_matrix[,6]==k,]))
    total_dist <- mean(dist(pc_matrix))
    ratio <- within_dist/total_dist
    out <- append(out,ratio)
  }
  return(1-mean(out))
}

glPca.plot<-function(pca_dat,pop_list,PrinX,PrinY){
  #PrinX<-as.integer(PrinX)
  #PrinY<-as.integer(PrinY)
  pca.scores<-as.data.frame(pca_dat$x)
  pca.scores$pop <- pop_list
  cols <- brewer.pal(n = length(unique(pop_list)), name = "BrBG")
  df<-as.data.frame(cbind(pca.scores[,PrinX],pca.scores[,PrinY]))
  df$Population<-pca.scores$pop
  colnames(df)<-c("PCX","PCY","Population")
  p <- ggplot(df, aes(x=df$PCX, y=PCY, colour="black",fill=Population))
  p <- p + geom_point(size=2,shape=21)
  var_frac <- pca_dat$sdev/sum(pca_dat$sdev^2)
  x<-signif(sum(dplyr::nth(var_frac,PrinX)) * 100, 3)
  y<-signif(sum(dplyr::nth(var_frac,PrinY)) * 100, 3)
  p <- p + xlab(paste(paste0("PC",PrinX),paste0("(",x,"%"), "of variance explained)",sep=" "))
  p <- p + ylab(paste(paste0("PC",PrinY),paste0("(",y,"%"), "of variance explained)",sep=" "))
  p <- p + geom_hline(yintercept = 0) 
  p <- p + geom_vline(xintercept = 0) 
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_text(size = 12))
  p <- p + theme(axis.title = element_text(size=14))
  p <- p + theme(legend.text = element_text(size = 12))
  p <- p + theme(legend.title = element_text(size = 14))
  #p <- p + theme(plot.margin = margin(2,.8,2,.8,"cm"))
  p  + scale_fill_manual(labels = unique(pop_list),values=brewer.pal(9,"Set1")) +
    scale_color_manual(values = "black")
}
library(RColorBrewer)
library(ggplot2)
return_miss_loc <- function(vcf, filt){
  dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
  q<-((apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) })) / rowSums(dp,na.rm = TRUE))
  q<-names(q[q >= filt])
  q
}

