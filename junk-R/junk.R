##' Count mismatches between all pairs of individuals
##'
##' @title mismatch.count.all
##' @inheritParams dups
##' @return a matrix, nrow(x) x nrow(y) with each entry the number of
##' mismatches if <tol, or tol
##' @author Chris Wallace
##' @export
mismatch.count.all <- function(x,y=NULL,tol=ncol(x)/100) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")
  
  counts <- matrix(as.integer(0),nrow=nrow(y),ncol=nrow(x),dimnames=list(rownames(y),rownames(x)))
  
  ret <- .C("samplediffall_c", x@.Data, y@.Data, counts=counts, as.integer(tol),
            nrow(x@.Data), ncol(x@.Data),
            nrow(y@.Data), ncol(y@.Data),
            PACKAGE="annotSnpStats")
  ## .Call("countmatches", Rx=x@.Data, Ry=y@.Data,
  ##   PACKAGE="annotSnpStats"))
  return(ret[["counts"]])
  
}
##' Count mismatches between specified pairs of individuals
##'
##' @title mismatch.count
##' @inheritParams dups
##' @return a vector with each entry the proportion of
##' mismatches if <tol, or tol
##' @examples
##' x<-example.data(1:10,1:10)
##' y<-x
##' y@@.Data[6,] <- as.raw("02")
##' mismatch.count(x,y)
##' @author Chris Wallace
mismatch.count <- function(x,y=NULL,tol=ncol(x)/100) {
  if(!identical(colnames(x), colnames(y)))
    stop("x and y need identical snps or sample comparison will be meaningless")
   if(!identical(rownames(x), rownames(y)))
    stop("x and y need identical samples for pairwise comparison.  Did you mean mismatch.count.all?")
  
  counts <- numeric(nrow(x))
  names(counts) <- rownames(x)
  
  ret <- .C("samplediff_c", x@.Data, y@.Data, counts=counts, as.integer(tol),
            nrow(x@.Data), ncol(x@.Data),
            PACKAGE="annotSnpStats")
  ## .Call("countmatches", Rx=x@.Data, Ry=y@.Data,
  ##   PACKAGE="annotSnpStats"))
  return(ret[["counts"]])
  
}
