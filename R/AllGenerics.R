##' get/set allele variables in an annotSnpStats object
##'
##' @title alleles
##' @param x object of class annotSnpStats
##' @param value character vector of length 2, giving the variable names for alelles 1 and 2
##' @return no return value
##' @author Chris Wallace
##' @export
##' @rdname alleles-methods
##' @keywords methods
setGeneric("alleles",
           def=function(x) {
             if(length(x@alleles))
               return(x@alleles)
             allele.options <- list(c("allele.1","allele.2"),
                                    c("allele.A","allele.B"),
                                    c("A1","A2"),
                                    c("a1","a2"))
             for(i in 1:length(allele.options)) {
               if(all(allele.options[[i]] %in% colnames(x@snps)))
                 return(allele.options[[i]])
             }
             warning("couldn't guess allele column names")
             return(NULL)
           })
##' @export
##' @rdname alleles-methods
##' @keywords methods
##' @aliases alleles<-,aSnpMatrix-method
setGeneric("alleles<-",
           def=function(x, value) {
             x@alleles <- value
             return(x)
           })


##' get phenotype variable in an annotSnpStats object
##'
##' @return character string 
##' @author Chris Wallace
##' @export
##' @rdname phenotype-methods
setGeneric("phenotype",
           def=function(x) {
             x@phenotype
           })
##' set phenotype variable in an annotSnpStats object
##'
##' @param x object of class aSnpMatrix or aXSnpMatrix
##' @param value character string naming phenotype column in samples(x)
##' @return no return value
##' @author Chris Wallace
##' @export
##' @rdname phenotype-methods
##' @aliases phenotype<-,aSnpMatrix,ANY-method
setGeneric("phenotype<-",
           def=function(x, value) {
             if(!(value %in% colnames(samples(x))))
               stop("phenotype must be name of a column in samples(x)")
             x@phenotype <- value
             return(x)
           })


##' Extract snps object from an object of class asnpMatrix
##'
##' @inheritParams sm
##' @export
##' @return snp information 
setGeneric("snps",
           def=function(obj) {
             return(obj@snps)
           })
##' Get/set samples object from an object of class asnpMatrix
##'
##' @inheritParams sm
##' @rdname samples-methods
##' @export
##' @return snp information 
setGeneric("samples",
           def=function(obj) {
             return(obj@samples)
           })
##' @export
##' @rdname samples-methods
##' @aliases samples<-,aSnpMatrix-method
setGeneric("samples<-",
           def=function(x, value) {
             if(nrow(value)!=nrow(x))
               stop("samples must have same nrow as x")
             if(!identical(rownames(value),rownames(x)))
               stop("samples must have identical rownames to x")
             x@samples <- value
             return(x)
           })


##' Add empty entries to x for SNPs found in y
##'
##' @rdname add-methods
##' @name addMethods
##' @param x aSnpMatrix or aXSnpMatrix
##' @param y aSnpMatrix or aXSnpMatrix
##' @export
##' @return aSnpMatrix: x with some empty (as.raw("00")) columns appended corresponding to SNPs in y
setGeneric("add.snps",function(x,y) standardGeneric("add.snps"))

##' Add empty entries to x for samples found in y
##'
##' @rdname add-methods
##' @name addMethods
##' @export
##' @return aSnpMatrix: x with some empty (as.raw("00")) rows appended corresponding to samples in y
setGeneric("add.samples",function(x,y) standardGeneric("add.samples"))

##' set rownames and colnames simulateously
##'
##' @export
##' @param x aSnpMatrix or aXSnpMatrix object
##' @param value list of two character vectors giving row and sample names
setGeneric("dimnames<-")
