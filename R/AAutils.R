##' Generate example annotSnpStats data objects
##'
##' Useful for testing the functions in annotSnpStats.  Note that the
##' example data are taken from the for.exercise dataset in snpStats
##' and contain 1000 samples, and 29501 SNPs, so don't expect to
##' get anything larger than that!
##'
##' @title example.data
##' @param m vector of sample indices
##' @param n vector of SNP indices
##' @return example data with given indices
##' @export
##' @examples
##' # 10 samples, 20 snps
##' example.data(1:10,1:20)
##' # different 10 samples, same 20 snps
##' example.data(101:110,1:20)
##' @author Chris Wallace
example.data <- function(m=1:100,n=1:10) {
  snps.10 <- snp.support <- subject.support <- NULL
  data("for.exercise", package="snpStats", envir=environment())
  if(length(m)==1)
    m <- 1:m
  if(length(n)==1)
    n <- 1:n
  if(any(m>nrow(snps.10)))
    stop("can only supply ",nrow(snps.10)," samples")
  if(any(n>nrow(snp.support)))
    stop("can only supply ",nrow(snp.support)," SNPs")
  snp.support$A1 <- as.character(snp.support$A1)
  snp.support$A2 <- as.character(snp.support$A2)
  new("aSnpMatrix",
      .Data=snps.10[m,n],
      snps=snp.support[n,,drop=FALSE],
      samples=subject.support[m,],
      phenotype="cc")
}
