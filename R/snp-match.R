snp.id <- function(x,col.x) {
  if(col.x==0) {
    return(rownames(x@snps))
  } else {
    return(x@snps[,col.x])
  }
}
sample.id <- function(x,col.x) {
  if(col.x==0) {
    return(rownames(x@samples))
  } else {
    return(x@samples[,col.x])
  }
}

##' Match snps in two aSnpMatrix objects
##'
##' Extract indices of overlapping SNPs in two aSnpMatrix objects.  By
##' default, match is on \code{rownames(x@@snps)}, \code{rownames(y@@snps)}.
##' @title snp.match
##' @param x aSnpMatrix object
##' @param y aSnpMatrix object
##' @param col.x optional, column identifier for \code{x@@snps} for matching
##' @param col.y optional, column identifier for \code{y@@snps} for matching
##' @return a named list of length two with elements
##' \describe{
##'
##' \item{x}{indices of matching SNPs in x}
##'
##' \item{y}{indices of matching SNPs in y}
##'
##' }
##' @author Chris Wallace
##' @export
snp.match <- function(x,y,col.x=0,col.y=0) {
  id.x <- snp.id(x,col.x)
  id.y <- snp.id(y,col.y)
  ids <- intersect(id.x,id.y)
  if(!length(ids))
    stop("no overlapping SNPs found")
  m.x <- match(ids,id.x)
  m.y <- match(ids,id.y)
  return(list(x=m.x, y=m.y))
}

##' Match samples in two aSnpMatrix objects
##'
##' Extract indices of overlapping samples in two aSnpMatrix objects.  By
##' default, match is on \code{rownames(x@@samples)}, \code{rownames(y@@snps)}.
##' @title snp.match
##' @param x aSnpMatrix object
##' @param y aSnpMatrix object
##' @param col.x optional, column identifier for \code{x@@samples} for matching
##' @param col.y optional, column identifier for \code{y@@samples} for matching
##' @return  a named list of length two with elements
##' \describe{
##'
##' \item{x}{indices of matching samples in x}
##'
##' \item{y}{indices of matching samples in y}
##'
##' }
##' @author Chris Wallace
##' @export
sample.match <- function(x,y,col.x=0,col.y=0) {
  id.x <- sample.id(x,col.x)
  id.y <- sample.id(y,col.y)
  ids <- intersect(id.x,id.y)
  if(!length(ids))
    stop("no overlapping samples found")
  m.x <- match(ids,id.x)
  m.y <- match(ids,id.y)
  return(list(x=m.x, y=m.y))
}
##' Count mismatches between pairs of individuals
##'
##' @title mismatch.count
##' @inheritParams dups
##' @return a matrix, nrow(x) x nrow(y) with each entry the number of
##' mismatches if <tol, or tol
##' @author Chris Wallace
##' @export
mismatch.count <- function(x,y=NULL,tol=ncol(x)/100) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")
  
  counts <- matrix(as.integer(0),nrow=nrow(y),ncol=nrow(x),dimnames=list(rownames(y),rownames(x)))
  
  ret <- .C("samplediff_c", x@.Data, y@.Data, counts=counts, as.integer(tol),
            nrow(x@.Data), ncol(x@.Data),
            nrow(y@.Data), ncol(y@.Data),
            PACKAGE="annotSnpStats")
  ## .Call("countmatches", Rx=x@.Data, Ry=y@.Data,
  ##   PACKAGE="annotSnpStats"))
  return(ret[["counts"]])
  
}
##' Find indices of possible sample duplications between two aSnpStats objects
##'
##' Each pair of samples from x and y are compared in turn.  If the
##' number of mismatched and non-missing genotypes exceeds tol, the
##' pair are assumed to be non-duplicates, and counting proceeds to
##' the next pair.  If the total number of mismatched and non-missing
##' genotypes is <tol, then the indices of the sample pair are stored,
##' and returned together with the number of mismatches and the number
##' of non-missing genotypes compared.
##' @title Find duplicate samples
##' @param x aSnpStats object
##' @param y aSnpStats object
##' @param tol maximum number of mismatched genotypes allowed for
##' duplicate samples
##' @param type by default, dups compares only homs vs hets, to allow
##' for differently labelled alleles.  Set type="all" to allow the two
##' kinds of homozygote genotypes to count as a mismatch
##' @return a matrix, with four columns: index of dup in x, index of
##' dup in y, number of mismatches, number of comparisons
##' @author Chris Wallace
##' @export
dups <- function(x,y,tol=ncol(x)/50,type=c("hethom","all")) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")
  type <- switch(match.arg(type),
                 all=0,
                 hethom=1)
  
  pBar <- txtProgressBar( min = 0, max = nrow(x) - 1, style = 3 )
  ret <- .Call("countdiffs", x@.Data, y@.Data, as.integer(min(tol,1)), as.integer(type), pBar,               
               PACKAGE="annotSnpStats")
  cat("\n") # end progress bar
  colnames(ret) <- c("index.x","index.y","mismatch","total")
  return(ret)
}
  
##' Drop SNPs from an annotSnpStats object according to SNP qc summary stats
##'
##' \code{thr} must be a named character vector, with names
##' corresponding to column names in the matrix returned by
##' \code{\link[snpStats]{col.summary}}.  The elements of the vector
##' should be character objects.  The first character must be '<' or
##' '>' and the remainder, on conversion to as.numeric, gives the
##' threshold for keeping SNPs.  \code{z.HWE} is special, and applies
##' to abs(z.HWE).
##' 
##' @title snp.trim
##' @param x annotSnpStats object
##' @param thr named character vector, see Details.
##' @examples
##' data(for.exercise, package="snpStats")
##' ## find rare SNPs with high call rates that do not deviate from HWE
##' snps.ok <- snp.trim(snps.10, thr=c(Call.rate=">0.99", MAF="<0.03", z.HWE="<3"))
##' @return an annotSnpStats object with only SNPs that meet the criteria specified.
##' @author Chris Wallace
##' @export
snp.trim <- function(x, thr=c(Call.rate=">0.99", MAF=">0.03", z.HWE="<5")) {
  cs <- col.summary(x)
  cs[,"z.HWE"] <- abs(cs[,"z.HWE"])
  cn <- colnames(cs)
  snps.ok <- TRUE
  for(i in 1:length(thr)) {
    col <- which(cn==names(thr)[i])
    if(!length(col)) {
      warning(paste("column",names(thr)[i],"not found in col.summary object"))
      next
    }
    direction=substr(thr[i],1,1)
    if(!(direction %in% c(">","<")))
      stop("direction must be '>' or '<'")
    threshold <- as.numeric(substr(thr[i],2,nchar(thr[i])))
    pass.this <- eval(call(direction,cs[,col],threshold))
    cat(names(thr)[i],direction,threshold,"\n")
    print(table(pass.this))
    snps.ok <- snps.ok & pass.this
  }
  cat("Overall\n")
  print(table(snps.ok))
  return(x[,which(snps.ok)])
}

complement <- function(str) {
  str <- tolower(str)
  str <- gsub("a","T",str)
  str <- gsub("t","A",str)
  str <- gsub("c","G",str)
  str <- gsub("g","C",str)
  return(str)
}

arev <- function(str) {
  ss <- strsplit(str,"/")
  ss <- lapply(ss, rev)
  sapply(ss, paste, collapse="/")  
}
##' switch alleles in x to match order in y
##'
##' @title align.alleles
##' @param x  annotSnpStats object
##' @param y  annotSnpStats object
##' @param do.plot
##' @export
##' @return  new annotSnpStats object derived from x, with alleles switched to match those in y
##' @author Chris Wallace
align.alleles <- function(x,y,do.plot=TRUE,mafdiff=0.01) {
  x.alleles <- apply(x@snps[,alleles(x)],1,paste,collapse="/")
  y.alleles <- apply(y@snps[,alleles(y)],1,paste,collapse="/")
  print(table(x.alleles, y.alleles))
  ok <- x.alleles == y.alleles | x.alleles == complement(y.alleles)
  sw <- x.alleles == arev(complement(y.alleles)) | x.alleles == arev(y.alleles)
  if(!all(ok | sw))
    stop("couldn't work out all switches")
  if(any(ok & sw)) {
   ind <- ok & sw
    cat(sum(ind),"SNPs have alleles not completely resolvable without strand information, confirming guess by checking allele freqs\n")

     x.cs <- col.summary(x[,ind])
    y.cs <- col.summary(y[,ind])

    rdiff <- x.cs[,"RAF"] - y.cs[,"RAF"]
    sw2 <- ifelse(abs(x.cs[,"RAF"] - y.cs[,"RAF"]) < abs(1 - x.cs[,"RAF"] - y.cs[,"RAF"]), FALSE, TRUE)
    sw2[ abs(x.cs[,"MAF"]-0.5)<mafdiff ] <- NA
    cat(sum(is.na(sw2)),"SNPs not resolvable (MAF too close to 0.5)\n")

    sw[ind] <- sw2
 }

  x@.Data[,is.na(sw)] <- as.raw("00")
  x <- switch.alleles(x, which(sw))
  x@snps[, alleles(x) ] <- y@snps[, alleles(y)]

  if(do.plot) {
    x.cs <- col.summary(x)
    y.cs <- col.summary(y)
    plot(x.cs[,"RAF"], y.cs[,"RAF"],main="RAF after switching",xlab="x",ylab="y")
  }
  
  return(x)
}


    ## plot(x.cs[,"RAF"],y.cs[,"RAF"])
    ## adiff <- paste(x.alleles[ind], y.alleles[ind])
    ## boxplot(rdiff ~ adiff)
    ## library(ggplot2)
    ## qplot(x.cs[,"RAF"],y.cs[,"RAF"],col=adiff, facets=~adiff)
    
    ## bpm <- read.table("/ipswich/data/Immunochip/Illumina-annotation-2010-09-07/Genotypes/Immuno_BeadChip_11419691_B.csv",sep=",",quote="",comment.char="",as.is=TRUE,header=TRUE)

    ## m <- match(colnames(y)[ind], make.names(bpm$Name))
    ## ilmn <- bpm[m,"IlmnID"]
    ## tr <- sub(".*([TB]_[RF]).*","\\1",ilmn)
    ## ss <- bpm[m,"SourceStrand"]
    ## is <- bpm[m,"IlmnStrand"]
    ##  qplot(x.cs[,"RAF"],y.cs[,"RAF"],col=adiff, facets=adiff~bpm[m,"IlmnStrand"])
    ##  qplot(x.cs[,"RAF"],y.cs[,"RAF"],col=adiff, facets=ss ~ is)
