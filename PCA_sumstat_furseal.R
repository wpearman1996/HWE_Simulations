library("optparse")
library(hierfstat)

option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="VCF name", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory to look in for VCF & store output", metavar="character"),
  make_option(c("-p", "--pre"), type="character", default=NULL, 
              help="output prefix name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$vcf)){
  print_help(opt_parser)
  stop("Input vcf is required", call.=FALSE)
}
if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory is required", call.=FALSE)
}
if (is.null(opt$pre)){
  print_help(opt_parser)
  stop("Output prefix is required", call.=FALSE)
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

library("adegenet")
library("pegas")
library("vcfR")
library("StAMPP")

setwd(opt$dir)
file <- read.vcfR(opt$vcf,verbose = FALSE)
file_gl <- vcfR2genlight(file)


pop_list<-read.table("./popmap.tsv")
file_gl$pop<-as.factor(pop_list$V2[match(file_gl$ind.names,pop_list$V1)])
pops <- as.character(unique(pop_list$V2))

file_gl_fst<-file_gl
fst_matrix<-stamppFst(file_gl_fst,nboots = 0)
file_gi <- gl2gi(file_gl)
file_fis<-basic.stats(genind2hierfstat(file_gi,pop=file_gl_fst$pop))
file_neis_dist<-genet.dist(genind2hierfstat(file_gi,pop=file_gl_fst$pop))
y<-scaleGen(file_gi,NA.method="mean",scale=F)
z<-prcomp(y,center=F,scale=F)$x
library(stringr)
vcfname<-word(opt$vcf,1,2,sep="[.]")
write.csv(file_fis$perloc,paste0(vcfname,"_summary_stats_perloc.csv"))
write.csv(mean(file_neis_dist),paste0(vcfname,"neis_dist.csv"))
write.csv(z,paste0(vcfname,"_PCA_matrix.csv"))
write.csv(fst_matrix,paste0(vcfname,"_Fst_matrix.csv"))
write.csv(file_fis$overall,paste0(vcfname,"_summary_stats.csv"))