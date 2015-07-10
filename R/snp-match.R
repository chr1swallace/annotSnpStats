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
##' default, match is on \code{rownames(snps(x))}, \code{rownames(snps(y))}.
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
##' @example examples/match.R
##' @author Chris Wallace
##' @export
snp.match <- function(x,y,col.x=0,col.y=0) {
  id.x <- snp.id(x,col.x)
  id.y <- snp.id(y,col.y)
  ids <- intersect(id.x,id.y)
  if(!length(ids)) {
    message("no overlapping SNPs found")
    return(NULL)
  }
  m.x <- match(ids,id.x)
  m.y <- match(ids,id.y)
  message(length(m.x)," matching SNPs found")
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
##' @example examples/match.R
##' @author Chris Wallace
##' @export
sample.match <- function(x,y,col.x=0,col.y=0) {
  id.x <- sample.id(x,col.x)
  id.y <- sample.id(y,col.y)
  ids <- intersect(id.x,id.y)
  if(!length(ids)) {
    message("no matching samples found")
    return(NULL)
  }
  message(length(ids)," matching samples found")
  m.x <- match(ids,id.x)
  m.y <- match(ids,id.y)
  return(list(x=m.x, y=m.y))
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
##' @param stopatone if TRUE, assume each sample in x can have at most
##' one match in y, and vice versa. This makes things faster, and
##' should be safe assuming x and y themselves contain no internal
##' duplicates so is set to TRUE by default, but set it to FALSE if
##' you want to catch multiple matches.
##' @return a matrix, with four columns: index of dup in x, index of
##' dup in y, number of mismatches, number of comparisons
##' @examples
##' ## example data where samples 6:10 in x are the same as 1:5 in y
##' x <- example.data(1:10,1:500)
##' y <- example.data(6:15,1:500)
##' dups(x,y)
##' @author Chris Wallace
##' @export
dups <- function(x,y,tol=ncol(x)/50,type=c("hethom","all"),stopatone=TRUE) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")
  type <- switch(match.arg(type),
                 all=0,
                 hethom=1)
  
#  pBar <- txtProgressBar( min = 0, max = nrow(x) - 1, style = 1 )
#  ret <- .Call("countdiffs", x@.Data, y@.Data, as.integer(max(tol,0)),
  ret <- .Call("annotSnpStats_dups", x@.Data, y@.Data, as.integer(max(tol,0)),
               as.integer(type), as.integer(stopatone), as.raw("00"), as.raw("02"),
               PACKAGE="annotSnpStats")
#  cat("\n") # end progress bar
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
##' @title snp.qc
##' @param x annotSnpStats object
##' @param thr named character vector, see Details.
##' @examples
##' X<-example.data(1000,1000)
##' X
##' ## find rare SNPs with high call rates that do not deviate from HWE
##' X.rare.qcpass <- snp.qc(X, thr=c(Call.rate=">0.99", MAF="<0.05", z.HWE="<3"))
##' X.rare.qcpass
##' @return an annotSnpStats object with only SNPs that meet the criteria specified.
##' @author Chris Wallace
##' @export
snp.qc <- function(x, thr=c(Call.rate=">0.99", MAF=">0.03", z.HWE="<5")) {
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

## complement <- function(str) {
##   str <- tolower(str)
##   str <- gsub("a","T",str)
##   str <- gsub("t","A",str)
##   str <- gsub("c","G",str)
##   str <- gsub("g","C",str)
##   return(str)
## }

## arev <- function(str) {
##   ss <- strsplit(str,"/")
##   ss <- lapply(ss, rev)
##   sapply(ss, paste, collapse="/")  
## }
##' Complement genotypes
##'
##' ie A/G -> T/C
##' 
##' @param x character vector of genotypes
##' @export
##' @return character vector of genotypes on the alternative strand
##' @examples
##' g.complement(c("A/G","A/T"))
g.complement <- function(x) {
  x <- toupper(x)
  switches <- c("A"="t","T"="a","C"="g","G"="c")
  for(i in seq_along(switches))
    x <- sub(names(switches)[i],switches[i],x)
  toupper(x)
}

##' Reverse alleles in a genotype
##'
##' ie A/G -> G/A
##'
##' @param x character vector of genotypes
##' @param sep character with which to separate alleles. Default is "/".
##' @export
##' @return character vector of reversed genotypes 
##' @examples
##' g.rev(c("A/G","A/T"))
g.rev <- function(x,sep="/") {
  sapply(strsplit(x,sep),function(g) paste(rev(g),collapse="/"))
}
##' count specific genotypes
##'
##' sum of numbers in tt indexed by cbind(xind,yind). Allows that not
##' all of xind and yind may be represented in tt
##'
##' @param tt table of genotype counts
##' @param xind row indices
##' @param yind y indices
##' @export
##' @return sum of numbers indexed by cbind(xind,yind)
g.count <- function(tt,xind,yind) {
  ind <- cbind(xind,yind)
  ind <- ind[xind %in% rownames(tt) & yind %in% colnames(tt), , drop=FALSE]
  sum(tt[ind],na.rm=TRUE)
}

count.switches <- function(tt) {
  genos <- c("A/C","A/G","C/A","C/T","G/A","G/T","T/C","T/G")
  rev.genos <- g.rev(genos)
  str.genos <- g.complement(genos)
  revstr.genos <- g.rev(str.genos)
  return(c(nochange=g.count(tt,genos,genos),
           rev = g.count(tt,genos,rev.genos),
           str.n=g.count(tt,genos,str.genos),
           revstr.n=g.count(tt,genos,revstr.genos)))
}
##' define possible allele switching classes
##'
##' @title g.class
##' @param x vector of allele codes from dataset X
##' @param y vector of allele codes from dataset Y, same length as x
##' @return character vector of allele switching classes
##' @export
##' @examples
##' alleles.X <- c(snp1="A/G",snp2="A/G",snp3="A/G",snp4="A/G",snp5="A/T",snp6="A/T")
##' alleles.Y <- c(snp1="A/G",snp2="G/A",snp3="T/C",snp4="C/T",snp5="A/T",snp6="T/A")
##' classes <- g.class(x=alleles.X,y=alleles.Y)
##' cbind(alleles.X,alleles.Y,classes)
##' @author Chris Wallace
g.class <- function(x,y) {
  if(!identical(names(x),names(y)))
    stop("x and y must relate to same SNPs")
  mat <- matrix(FALSE,length(x),4,dimnames=list(names(x),c("nochange","rev","comp","revcomp")))
  ## nochange
  mat[ , "nochange" ] <- x==y
  mat[, "rev"] <- x==g.rev(y)
  mat[,"comp"] <- x==g.complement(y)
  mat[,"revcomp"] <- x==g.rev(g.complement(y))
  indels <- x %in% c("I/D","D/I")
  if(any(indels))
    mat[indels,c("comp","revcomp")] <- FALSE
  ret <- character(nrow(mat))
  rs <- rowSums(mat)
  if(length(wh <- which(rs>1))) # ambiguity first
    ret[wh] <- "ambig"  
  if(length(wh <- which(rs==0))) # impossible
    ret[wh] <- "impossible"
  if(length(wh <- which(rs==1))) # impossible
    ret[wh] <- colnames(mat)[ apply(mat[wh,,drop=FALSE],1,which) ]
  return(ret)
}

##' switch alleles in x to match order in y
##'
##' @title align.alleles
##' @param x annotSnpStats object
##' @param y annotSnpStats object
##' @param do.plot if TRUE (default), generate a summary plot that can be a useful visual check nothing has gone wrong.  The points should lie close to the line of equality.
##' @param mafdiff SNPs with MAF within mafdiff of 0.5 will not be
##' aligned automatically.  This is of concern only for T/A and C/G
##' SNPs, which are not uniquely resolveable by allele codes.
##' @param known.dups if duplicate samples exist, supply a list
##' returned by running dups(x,y).  The alignment of genotypes in
##' these duplicate samples will be used to inform alignment of
##' genotypes in other samples.
##' @examples
##' x <- example.data(1:100,1:10)
##' y <- example.data(101:200,1:10)
##' ## switch columns 6:10
##' y.switched <- switch.alleles(y,6:10)
##' ## automatically switch back
##' y.aligned <- align.alleles(y.switched,x,do.plot=TRUE)
##' ## check by comparing counted allele frequencies for SNPs 1:5 and 6:10
##' cbind(x=col.summary(x)[,"RAF"],
##'       y=col.summary(y)[,"RAF"], 
##'       y.switched=col.summary(y.switched)[,"RAF"],
##'       y.aligned=col.summary(y.aligned)[,"RAF"])
##' @export
##' @return  new annotSnpStats object derived from x, with alleles switched to match those in y
##' @author Chris Wallace
align.alleles <- function(x,y,do.plot=TRUE,mafdiff=0.1,known.dups=NULL) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y should contain the same SNPs in the same order")
  if(any(is.na(x@snps[,alleles(x)])) || any(is.na(y@snps[,alleles(y)]))) {
    message("Missing alleles found, trying to fill in.  Missingness table before correction:")
    x.missing <- apply(is.na(x@snps[,alleles(x)]),1,any)
    y.missing <- apply(is.na(y@snps[,alleles(y)]),1,any)
    print(table(x.missing,y.missing))
    x1 <- x@snps[,alleles(x)[1]]
    x2 <- x@snps[,alleles(x)[2]]
    y1 <- y@snps[,alleles(y)[1]]
    y2 <- y@snps[,alleles(y)[2]]    
    m <- is.na(x1) & !is.na(x2) & !y.missing
    if(any(m)) {
      x1[ which(m & x2==y2) ] <- y1[ which(m & x2==y2) ]
      x1[ which(m & x2==y1) ] <- y2[ which(m & x2==y1) ]
    }    
    m <- !is.na(x1) & is.na(x2) & !y.missing
    if(any(m)) {
      x2[ which(m & x1==y2) ] <- y1[ which(m & x1==y2) ]
      x2[ which(m & x1==y1) ] <- y2[ which(m & x1==y1) ]
    }
    m <- is.na(y1) & !is.na(y2) & !x.missing
    if(any(m)) {
      y1[ which(m & y2==x2) ] <- x1[ which(m & y2==x2) ]
      y1[ which(m & y2==x1) ] <- x2[ which(m & y2==x1) ]
    }    
    m <- !is.na(y1) & is.na(y2) & !x.missing
    if(any(m)) {
      y2[ which(m & y1==x2) ] <- x1[ which(m & y1==x2) ]
      y2[ which(m & y1==x1) ] <- x2[ which(m & y1==x1) ]
    }
    x@snps[,alleles(x)[1]] <-     x1
    x@snps[,alleles(x)[2]] <-     x2
    y@snps[,alleles(y)[1]] <-     y1
    y@snps[,alleles(y)[2]] <-     y2
    x.missing <- apply(is.na(x@snps[,alleles(x)]),1,any)
    y.missing <- apply(is.na(y@snps[,alleles(y)]),1,any)
    message("Missingness table after correction:")
    print(table(x.missing,y.missing)) 
  }

  x.alleles <- apply(x@snps[,alleles(x)],1,paste,collapse="/")
  y.alleles <- apply(y@snps[,alleles(y)],1,paste,collapse="/")
  message("These are the allele codes as currently defined, before any switching:")
  print(tt <- as.matrix(table(x.alleles, y.alleles)))
  ## genotype classes
  sw.class <- g.class(x.alleles,y.alleles)
  any.comp <- any(sw.class %in% c("comp","revcomp"))
  any.ambig <- any(sw.class=="ambig")
  sw <- sw.class %in% c("rev","revcomp")
  if(length(wh <- which(sw.class=="impossible"))) {
    message(length(wh)," impossible genotype calls found. These genotypes will be set to missing.")
    sw[wh] <- NA
  }
  sw[sw.class=="impossible"] <- NA
  if(any.comp & any.ambig) { # there are reverse complements in the distinguishable cases
    ind <- which(sw.class=="ambig")
    message(length(ind)," SNPs have alleles not completely resolvable without strand information, confirming guess by checking allele freqs.")
    x.cs <- col.summary(x[,ind])
    y.cs <- col.summary(y[,ind])

    rdiff <- x.cs[,"RAF"] - y.cs[,"RAF"]
    sw2 <- ifelse(abs(x.cs[,"RAF"] - y.cs[,"RAF"]) < abs(1 - x.cs[,"RAF"] - y.cs[,"RAF"]), FALSE, TRUE)
    too.close <- abs(x.cs[,"MAF"]-0.5)<mafdiff
    if(any(too.close) & !is.null(known.dups)) {
      message(sum(too.close), "/", length(ind)," ambiguous SNPs close to 50% MAF.  Using known dups to resolve")
      if(is.character(known.dups[,1])) {
        m1 <- match(known.dups[,1],rownames(x))
        m2 <- match(known.dups[,2],rownames(y))
      } else {
        m1 <- known.dups[,1]
        m2 <- known.dups[,2]        
      }
      if(length(wh <- which(duplicated(m1)))) {
        m1 <- m1[-wh]
        m2 <- m2[-wh]
      } 
      if(length(wh <- which(duplicated(m2)))) {
        m1 <- m1[-wh]
        m2 <- m2[-wh]
      } 
      
      xn <- as(x[m1,ind[too.close]], "numeric")
      yn <- as(y[m2,ind[too.close]], "numeric")
      cor.asis <- sapply(seq_along(which(too.close)), function(i)
                         cor(xn[,i], yn[,i], use="pair"))
      
      cor.asis[ cor.asis > -0.8 & cor.asis < 0.8 ] <- NA # NA unless correlation is pretty strong
      sw2[too.close] <- cor.asis < 0
      too.close <- too.close[is.na(cor.asis)]
    }
    
    if(any(too.close)) {
      can.match <- sw.class %in% c("comp","nochange","rev","revcomp")
      xsw <- switch.alleles(x[,-ind],sw[-ind])
      ysw <- y[,-ind]
      ## step through too.close SNPs checking signed correlation
      message("using signed correlation for ",sum(too.close)," SNPs too close to 50% MAF")
      ldx <- ld(xsw,x[,ind[too.close],drop=FALSE], stats="R")
      ldy <- ld(ysw,y[,ind[too.close],drop=FALSE], stats="R")
      ldx[abs(ldx)<0.04] <- NA ## drop uncorrelated - have no information
      ldy[abs(ldy)<0.04] <- NA ## drop uncorrelated - have no information
      cor.sw <- sapply(1:ncol(ldx), function(j) cor(ldx[,j], ldy[,j], use="pair"))
      cor.sw[ abs(cor.sw)<0.8 ] <- NA # NA unless correlation is pretty strong
      sw2[too.close] <- cor.sw < 0
      too.close <- too.close[is.na(cor.sw)]
    }
    message(sum(is.na(sw2))," SNPs not resolvable (MAF too close to 0.5).")
    sw[ind] <- sw2
  }

  if(!any.comp & any.ambig) { # there are no reverse complements in distinguishable cases
    ind <- which(sw.class=="ambig")
    message(length(ind)," SNPs have alleles not completely resolvable without strand information,\nbut there is no evidence of strand switches amongst SNPs which are resolvable.\nAssuming fixed strand.")
     ind <- which(sw.class=="ambig")
     sw2 <- x.alleles[ind]==g.rev(y.alleles[ind])
     sw[ind] <- sw2
  }

  if(any(is.na(sw)))
    x@.Data[,is.na(sw)] <- as.raw("00")
  if(length(wh <- which(sw))) {
    x <- switch.alleles(x, wh)
    x@snps[wh, alleles(x) ] <- y@snps[wh, alleles(y)]
  }

  if(do.plot) {
    x.cs <- col.summary(x)
    y.cs <- col.summary(y)
    plot(x.cs[,"RAF"], y.cs[,"RAF"],main="RAF after switching",xlab="x",ylab="y",pch="+")
    abline(0,1,col="red",lty=1)
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
