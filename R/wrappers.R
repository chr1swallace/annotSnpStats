## wrappers for S3 functions in snpMatrix that don't obey inheritances
##' Wrapper to allow snpStats::ld to be called with aSnpMatrix objects.
##'
##' See \code{\link[snpStats]{ld}} for help on the function itself.
##'
##' @export
##' @param x object of class aSnpMatrix or aXSnpMatrix
##' @param y object of same class as x
##' @param ... other arguments passed to \code{snpStats::ld()}
ld <- function(x, y=NULL, ...) {
  if(is(x,"aSnpMatrix"))
    x <- sm(x)
  if(!is.null(y) && is(y,"aSnpMatrix"))
    y <- sm(y)
  snpStats::ld(x=x, y=y, ...)
}

