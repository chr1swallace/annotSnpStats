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

dups <- function(x,y=NULL) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")

  counts <- matrix(as.integer(0),nrow=nrow(y),ncol=nrow(x),dimnames=list(rownames(y),rownames(x)))
  
  ret <- .C("samplecmp", x@.Data, y@.Data, counts=counts,
            nrow(x@.Data), ncol(x@.Data),
            nrow(y@.Data), ncol(y@.Data),
            PACKAGE="annotSnpStats")
  return(ret[["counts"]])
  
}

dups.c <- function(x,y=NULL,tol=ncol(x)/50,type=c("diff","same")) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")

  counts <- matrix(as.integer(0),nrow=nrow(y),ncol=nrow(x),dimnames=list(rownames(y),rownames(x)))

  type <- match.arg(type)
 
  ret <- switch(type,
                ## "diff" = .Call("countdiffs", Rx=x@.Data, Ry=y@.Data, maxDiff=maxdiff,
                ##   PACKAGE="annotSnpStats"),
                "diff" =.C("samplediff_c", x@.Data, y@.Data, counts=counts, as.integer(tol),
                  nrow(x@.Data), ncol(x@.Data),
                  nrow(y@.Data), ncol(y@.Data),
                  PACKAGE="annotSnpStats"),
                "same" = .C("samplecmp", x@.Data, y@.Data, counts=counts,
                  nrow(x@.Data), ncol(x@.Data),
                  nrow(y@.Data), ncol(y@.Data),
                  PACKAGE="annotSnpStats"))
                ## .Call("countmatches", Rx=x@.Data, Ry=y@.Data,
                ##   PACKAGE="annotSnpStats"))
  return(ret[["counts"]])
  
}

dups <- function(x,y,tol=ncol(x)/50) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")

  ret <- .Call("countdiffs", x@.Data, y@.Data, as.integer(tol),
                  PACKAGE="annotSnpStats")
  colnames(ret) <- c("index.x","index.y","mismatch","total")
  return(ret)
  
}
