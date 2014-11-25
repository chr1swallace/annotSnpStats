#' Annotated SnpMatrix class
#'
#' SnpMatrix or XSnpMatrix objects, tied to sample and SNP support objects
#'
#' This is a small class to allow subsetting and binding operations to be applied to SnpMatrix objects together with their SNP and sample support objects in single operations.  It means the different objects are assured to always line up.
#'
#'@slot .Data \code{SnpMatrix} or \code{XSnpMatrix} object
#'@slot snps \code{data.frame} with rownames matching colnames of \code{.Data}
#'@slot samples \code{data.frame} with rownames matching rownames of \code{.Data}
#'@slot diploid logical vector of length == \code{ncol{.Data}}, TRUE indicates sample from a diploid individual
#'
#' @exportClass aSnpMatrix
#' @author Chris Wallace
setClass("aSnpMatrix",
         representation(snps="data.frame",samples="data.frame",phenotype="character",alleles="character"),
         contains="SnpMatrix",
         validity=function(object) {
           if(ncol(object@.Data) && !identical(rownames(object@snps),colnames(object@.Data)))
             stop("rownames of snps object must exactly match colnames of SnpMatrix object")
           if(nrow(object@.Data) && !identical(rownames(object@samples),rownames(object@.Data)))
             stop("rownames of samples object must exactly match rownames of SnpMatrix object")
           if(length(object@phenotype) && !(object@phenotype %in% colnames(object@samples)))
             stop("phenotype must specify a column name of object@samples")
           if(length(object@alleles) && (!all(object@alleles %in% colnames(object@snps)) ||
                                         length(object@alleles)!=2))
             stop("alleles must specify a two columns of object@snps")
         })

#' @rdname aSnpMatrix-class
#' @exportClass aXSnpMatrix
setClass("aXSnpMatrix",
         representation(snps="data.frame",samples="data.frame",phenotype="character",alleles="character",dipload="logical"),
         ## should use contains=XSnpMatrix, but bug in snpStats prevents this
         contains="SnpMatrix",
         validity=function(object) {
           if(ncol(object@.Data) && !identical(rownames(object@snps),colnames(object@.Data)))
             stop("rownames of snps object must exactly match colnames of SnpMatrix object")
           if(nrow(object@.Data) && !identical(rownames(object@samples),rownames(object@.Data)))
             stop("rownames of samples object must exactly match rownames of SnpMatrix object")
           if(length(object@phenotype) && !(object@phenotype %in% colnames(object@samples)))
             stop("phenotype must specify a column name of object@samples")
           if(length(object@alleles) && (
                      !all(object@alleles %in% colnames(object@snps)) ||
                      length(object@alleles)!=2))
             stop("alleles must specify a two columns of object@snps")
         })
