setwd("/home/peawi142/HWE_Simulations/")
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
library(ggplot2)
library(ggridges)
library(StAMPP)
library(stringr)
library(pegas)
library(vcfR)
library(dartR)
library(RColorBrewer)
library(pophelper)
admix_plot<-function(clump_align_tab,nrep,nind,k,axis,col_pal,name){
  line_locs<-c(30,60,90,120,150)
  z<-clump_align_tab
  barplot(t(z),col=col_pal,xlab="Admixture Proportions",ylab=NULL,
          space= 0,border=NA,axisnames = axis,horiz=T, main = name,
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
run_structure_analysis <- function(directory, k, pop_list,simulation,useclumpp){
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
  read_structure <- function(file, nind,k){
    x<-read.delim(file,
                  skip=(49+(k*3)),head=F)
    x<-x$V1[1:nind]
    library(stringr)
    x<-gsub("\\s+", " ", str_trim(x))
    x <- data.frame(do.call(rbind, strsplit(x, " ", fixed=TRUE)))
    x<-x[c(5:(4+k))]
  }
  filt<-c("nohwe","out_across","out_all","out_any")
  library("pophelper")
  for(i in 1:4){
    files<-paste0(directory,
                  list.files(pattern="*_f",path = directory))
    files<-files[grepl(paste0("*",filt[i],"*"),files)]
    files<-files[grepl(paste0("*K",k),files)]
    slist <-readQ(files=files)
    files<-slist
    if(useclumpp == TRUE){
      
      clumppExport(slist,useexe = useclumpp,dir_name = paste0("K",k,"_",filt[i],simulation))}}
  col_pal<-brewer.pal(k,"Paired")
  files<-vector()
  for ( i in 1:4){
    files[i]<-list.files(pattern = "*merged.txt",path=paste0("K",k,"_",filt[i],simulation))
    files[i]<-paste0(paste0("K",k,"_",filt[i],simulation,"/",files[i]))
  }
  
  slist <-readQ(files=files)
  names(slist)<-paste0(simulation,filt)
  clumppExport(slist,useexe = useclumpp,dir_name = "clumpped_filtered",)
  clumpped_dat<-read.table("./clumpped_filtered/pop_K6-combined-aligned.txt")
  clumpped_dat<-clumpped_dat[,2:7]
  clumpped_dat_sep <- list()
  clumpped_dat_sep$nofilt<-clumpped_dat[1:180,]
  clumpped_dat_sep$out_across<-clumpped_dat[181:360,]
  clumpped_dat_sep$out_all<-clumpped_dat[361:540,]
  clumpped_dat_sep$out_any<-clumpped_dat[541:720,]
  clumpped_dat_sep
}
gl2gi <- function(gl, probar=FALSE, verbose=NULL) {
  
  
  # FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(gl@other$verbose)) verbose=gl@other$verbose
  if (is.null(verbose)) verbose=2
  
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose>0) {
    cat("Start conversion....\n")
    ptm <- proc.time()[3]
    cat("Please note conversion of bigger data sets will take some time!\n" )
    cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
  }
  #convert to genind....
  x <- as.matrix(gl[,])
  
  if (probar) {pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)}
  
  if (is.null(gl@loc.all))  {
    gl@loc.all <- rep("A/T", nLoc(gl))
    gl@loc.all[1]<- "C/G"
  }
  
  
  homs1 <- paste(substr(gl@loc.all,1,1),"/",substr(gl@loc.all,1,1), sep = "")
  hets <-  gl@loc.all
  homs2 <- paste(substr(gl@loc.all,3,3),"/",substr(gl@loc.all,3,3), sep = "")
  xx <- matrix(NA, ncol=ncol(x), nrow=nrow(x))
  for (i in 1:nrow(x))
  {
    for (ii in 1:ncol(x))
    {
      
      inp <- x[i,ii]
      if (!is.na(inp))
      {
        
        if (inp==0) xx[i,ii] <- homs1[ii] else if (inp==1) xx[i,ii] <- hets[ii] else if (inp==2) xx[i,ii] <- homs2[ii]
      } else xx[i,ii]="-/-"
    }
    if (probar) {setTxtProgressBar(pb, i/nrow(x))}
    
    if (verbose==1) {
      cat("\nMatrix converted.. Prepare genind object...\n")}
  }
  if (probar) {close(pb)}
  
  gen<-df2genind(xx[,], sep="/", ncode=1, ind.names=gl@ind.names, pop = gl@pop, ploidy=2,  NA.char = "-")#, probar=probar)
  gen@other <- gl@other
  locNames(gen)<- locNames(gl)
  
  if (verbose==1)cat(paste("Finished! Took", round(proc.time()[3]-ptm),"seconds.\n") )
  
  return(gen)
  
}
gl.drop.pop <- function(x, pop.list, as.pop=NULL, recalc=FALSE, mono.rm=FALSE, verbose=NULL){
  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
  # SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # Population labels assigned?
  if(is.null(as.pop)){
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      cat("Fatal Error: Population assignments not detected, running compliance check\n")
      x <- gl.compliance.check(x,verbose=0)
    }
  }
  
  # Assign the new population list if as.pop is specified
  pop.hold <- pop(x)
  if (!is.null(as.pop)){    
    if(as.pop %in% names(x@other$ind.metrics)){
      pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
      if (verbose >= 2) {cat("  Temporarily setting population assignments to",as.pop,"as specified by the as.pop parameter\n")}
    } else {
      stop("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n")
    }
  }
  
  if (verbose >= 2) {
    cat("  Checking for presence of nominated populations\n")
  }
  for (case in pop.list){
    if (!(case%in%popNames(x))){
      if(verbose >= 1){cat("  Warning: Listed population",case,"not present in the dataset -- ignored\n")}
      pop.list <- pop.list[!(pop.list==case)]
    }
  }
  if (length(pop.list) == 0) {
    stop("Fatal Error: no populations listed to drop!\n")
  }
  
  # DO THE JOB
  
  # Remove populations
  
  if (verbose >= 2) {
    cat("  Deleting populations", pop.list, "\n")
  }
  
  # Delete listed populations, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
  x2 <- x[!x$pop%in%pop.list]
  pop.hold <- pop.hold[!x$pop%in%pop.list]
  x <- x2
  
  # Monomorphic loci may have been created
  x@other$loc.metrics.flags$monomorphs == FALSE
  
  # Remove monomorphic loci
  if(mono.rm){
    if(verbose >= 2){cat("  Deleting monomorphic loc\n")}
    x <- gl.filter.monomorphs(x,verbose=0)
  } 
  # Check monomorphs have been removed
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat("  Warning: Resultant dataset may contain monomorphic loci\n")
    }  
  }
  
  # Recalculate statistics
  if (recalc) {
    x <- gl.recalc.metrics(x,verbose=0)
    if(verbose >= 2){cat("  Recalculating locus metrics\n")}
  } else {
    if(verbose >= 2){
      cat("  Locus metrics not recalculated\n")
      x <- utils.reset.flags(x,verbose=0)
    }
  }
  
  # REPORT A SUMMARY
  
  if (verbose >= 3) {
    if (!is.null(as.pop)) {
      cat("  Summary of recoded dataset\n")
      cat(paste("    No. of loci:",nLoc(x),"\n"))
      cat(paste("    No. of individuals:", nInd(x),"\n"))
      cat(paste("    No. of levels of",as.pop,"remaining: ",nPop(x),"\n"))
      cat(paste("    No. of populations: ",length(unique((pop.hold))),"\n"))
    } else {
      cat("  Summary of recoded dataset\n")
      cat(paste("    No. of loci:",nLoc(x),"\n"))
      cat(paste("    No. of individuals:", nInd(x),"\n"))
      cat(paste("    No. of populations: ", nPop(x),"\n"))
    }  
  }
  
  # Reassign the initial population list if as.pop is specified
  
  if (!is.null(as.pop)){
    pop(x) <- pop.hold
    if (verbose >= 3) {cat("  Resetting population assignments to initial state\n")}
  }
  
  # ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call() 
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:", funname, "\n")
  }
  
  return(x)
}
gl.filter.monomorphs <- function (x, verbose=NULL) {
  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
  # SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
  # DO THE JOB
  
  # mml <- !( colMeans(as.matrix(x), na.rm=TRUE)%%2 == 0) ; Code for readability
  
  hold <- x
  na.counter <- 0
  loc.list <- array(NA,nLoc(x))
  
  if (verbose >= 2){
    cat("Identifying monomorphic loci\n")
  }  
  
  # Tag presence/absence data
  if (data.type=="SilicoDArT"){
    nL <- nLoc(x)
    matrix <- as.matrix(x)
    l.names <- locNames(x)
    for (i in 1:nL){
      row <- matrix[,i] # Row for each locus
      if (all(row == 0, na.rm=TRUE) | all(row == 1, na.rm=TRUE) | all(is.na(row))){
        loc.list[i] <- l.names[i]
        if (all(is.na(row))){
          na.counter = na.counter + 1
        }
      }
    }                          
  } 
  
  # SNP data
  if (data.type=="SNP"){
    nL <- nLoc(x)
    matrix <- as.matrix(x)
    lN <- locNames(x)
    for (i in 1:nL){
      row <- matrix[,i] # Row for each locus
      if (all(row == 0, na.rm=TRUE) | all(row == 2, na.rm=TRUE) | all(is.na(row))){
        loc.list[i] <- lN[i]
        if (all(is.na(row))){
          na.counter = na.counter + 1
        }
      }
    }                          
  } 
  
  # Remove NAs from list of monomorphic loci and loci with all NAs
  loc.list <- loc.list[!is.na(loc.list)]
  
  # remove monomorphic loc and loci with all NAs
  
  if(length(loc.list > 0)){
    if (verbose >= 2){    cat("  Removing monomorphic loci\n")} 
    x <- gl.drop.loc(x,loc.list=loc.list,verbose=0)
  } else {
    if (verbose >= 2){cat("  No monomorphic loci to remove\n")}
  }
  
  # Report results
  if (verbose >= 3) {
    cat("  Original No. of loci:",nLoc(hold),"\n")
    cat("  Monomorphic loci:", nLoc(hold)-nLoc(x)-na.counter,"\n")
    cat("  Loci scored all NA:",na.counter,"\n")
    cat("  No. of loci deleted:",nLoc(hold)-nLoc(x),"\n")
    cat("  No. of loci retained:",nLoc(x),"\n")
    cat("  No. of individuals:",nInd(x),"\n")
    cat("  No. of populations:",nPop(x),"\n")
  }
  
  # RESET THE FLAG
  
  x@other$loc.metrics.flags$monomorphs <- TRUE
  
  # ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose >= 1){
    cat("Completed:",funname,"\n")
  }  
  
  return (x)
  
}
gl.drop.loc <- function(x, loc.list=NULL, first=NULL, last=NULL, verbose=NULL){
  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  hold <- x
  
  # SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 1){
    if (all(x@ploidy == 1)){
      cat("  Processing Presence/Absence (SilicoDArT) data\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (!is.null(loc.list) && !is.null(first)){
    flag <- 'both'
    if (verbose >= 2){
      cat("  Both a range of loci and a list of loci to keep has been specified\n")
    } 
  } else if (!is.null(loc.list)){
    flag <- 'list'
    if (verbose >= 2){
      cat("  List of loci to keep has been specified\n")
    } 
  } else if (!is.null(first)){
    flag <- 'range'
    if (verbose >= 2){
      cat("  Range of loci to keep has been specified\n")
    } 
  } else {
    stop("  Fatal Error: Need to specify either a range of loci to keep, or specific loci to keep\n")
  }
  
  if (flag=='both' || flag=='list'){
    for (case in loc.list){
      if (!(case%in%locNames(x))){
        cat("  Warning: Listed loci",case,"not present in the dataset -- ignored\n")
        loc.list <- loc.list[!(loc.list==case)]
      }
    }
  }
  
  if (flag=='range'){
    if (first <=0){
      cat("  Warning: Lower limit to range of loci cannot be less than 1, set to 1\n)")
      first <- 1
    }
    if (first > nLoc(x)){
      cat("  Warning: Upper limit to range of loci cannot be greater than the number of loci, set to",nLoc(x),"\n)")
      last <- nLoc(x)
    }
    if (first > last){
      cat("  Warning: Upper limit is smaller than lower limit, reversed\n")
      tmp <- first
      first <- last
      last <- tmp
    }
  }
  
  # DO THE JOB
  
  if (verbose >= 2) {
    cat("    Deleteing the specified loci\n")
  }
  
  # Remove duplicated loci if specified
  
  if (!is.null(first) && !is.null(loc.list)){
    list.from.range <- locNames(x)[first:last]
    loc.list <- unique(c(loc.list,list.from.range))
  } else if (!is.null(first)) {
    loc.list <- locNames(x)[first:last]
  }
  if (length(loc.list) == 0) {
    cat("  Warning: no loci listed to delete! Genlight object returned unchanged\n")
    x2 <- x
  } else {
    # Remove loci flagged for deletion
    x2 <- x[,!x$loc.names%in%loc.list]
    x2@other$loc.metrics <- x@other$loc.metrics[!x$loc.names%in%loc.list,]
  }  
  
  # REPORT A SUMMARY
  
  if (verbose >= 3) {
    cat("\n  Summary of recoded dataset\n")
    cat(paste("    Original No. of loci:",nLoc(hold),"\n"))
    cat(paste("    No. of loci deleted:",nLoc(hold)-nLoc(x2),"\n"))
    cat(paste("    No. of loci retained:",nLoc(x2),"\n"))
    cat(paste("    No. of individuals:", nInd(x2),"\n"))
    cat(paste("    No. of populations: ", nPop(x2),"\n\n"))
  }
   
  # ADD TO HISTORY    
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose >= 1){  
    cat("Completed:",funname,"\n")
  }
  
  return(x2)
  
}
