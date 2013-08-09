##' Convert output of \code{read.plink()} to an \code{aSnpMatrix} object
##'
##' @title annot.plink
##' @param obj output of \code{read.plink()}
##' @return aSnpMatrix object
##' @author Chris Wallace
##' @export
annot.plink <- function(obj) {
  new("aSnpMatrix",
      .Data=obj$genotypes,
      snps=obj$map,
      samples=obj$fam)
}
