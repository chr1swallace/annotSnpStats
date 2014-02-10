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
      samples=obj$fam,
      alleles=c("allele.1","allele.2"),
      phenotype="affected")
}

##' Read plink input files and convert to an aSnpMatrix object
##'
##' @title annot.read.plink
##' @param filestub character object specifying the stub of a plink
##' file.  Files with .bed .bim, and .fam extensions will be read
##' using snpStats::read.plink and the result converted to an
##' aSnpMatrix object
##' @return aSnpMatrix object
##' @author Chris Wallace
##' @export
annot.read.plink <- function(filestub) {
  obj <- read.plink(fam=sprintf("%s.fam",filestub),
                       bed=sprintf("%s.bed",filestub),
                       bim=sprintf("%s.bim",filestub))
  new("aSnpMatrix",
      .Data=obj$genotypes,
      snps=obj$map,
      samples=obj$fam,
      alleles=c("allele.1","allele.2"),
      phenotype="affected")
}

setAs("aSnpMatrix", "SnpMatrix", function(from) new("SnpMatrix", from@.Data))
##' Extract snpMatrix object from an object of class aSnpMatrix
##'
##' @param obj object of class aSnpMatrix
##' @return object of class snpMatrix
##' @export
sm <- function(obj) {
  as(obj, "SnpMatrix")
}
##' Extract snps object from an object of class snpMatrix
##'
##' @inheritParams sm
##' @export
##' @return snp information 
snps <- function(obj) {
  obj@snps
}
##' Extract samples object from an object of class snpMatrix
##'
##' @inheritParams sm
##' @export
##' @return snp information 
samples <- function(obj) {
  obj@samples
}
