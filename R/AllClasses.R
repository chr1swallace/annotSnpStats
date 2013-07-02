#' Annotated SnpMatrix class
#'
#' SnpMatrix objects, tied to sample and SNP support objects
#'
#' This is a small class to allow subsetting and binding operations to be applied to SnpMatrix objects together with their SNP and sample support objects in single operations.  It means the different objects are assured to always line up.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{.Data}:}{\code{SnpMatrix} object}
#'    \item{\code{snps}:}{\code{data.frame} with rownames matching colnames of \code{.Data}}
#'    \item{\code{samples}:}{\code{data.frame} with rownames matching rownames of \code{.Data}}
#'  }
#'
#' @name aSnpMatrix
#' @rdname aSnpMatrix
#' @aliases aSnpMatrix-class
#' @exportClass aSnpMatrix
#' @author Chris Wallace
setClass("aSnpMatrix",
         representation(snps="data.frame",samples="data.frame",phenotype="character",alleles="character"),
         contains="SnpMatrix",
         validity=function(object) {
           if(!identical(rownames(object@snps),colnames(object@.Data)))
             stop("rownames of snps object must exactly match colnames of SnpMatrix object")
           if(!identical(rownames(object@samples),rownames(object@.Data)))
             stop("rownames of samples object must exactly match rownames of SnpMatrix object")
           if(length(object@phenotype) && !(object@phenotype %in% colnames(object@samples)))
             stop("phenotype must specify a column name of object@samples")
           if(length(object@alleles) && (
                      !all(object@alleles %in% colnames(object@snps)) ||
                      length(object@alleles)!=2))
             stop("alleles must specify a two columns of object@snps")
         })
