##' get/set allele variables in an annotSnpStats object
##'
##' @title alleles
##' @param object object of class annotSnpStats
##' @return no return value
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname alleles-methods
##' @aliases alleles-methods
##' @keywords methods
setGeneric("alleles",
           def=function(x, ..., value) {
             if(length(x@alleles))
               return(x@alleles)
             allele.options <- list(c("allele.1","allele.2"),
                             c("allele.A","allele.B"))
             for(i in 1:length(allele.options)) {
               if(all(allele.options[[i]] %in% colnames(x@snps)))
                 return(allele.options[[i]])
             }
             warning("couldn't guess allele column names")
             return(NULL)
           })
##' @param ... character vector of length 2, giving the variable names for alelles 1 and 2
##' @export
##' @docType methods
##' @rdname alleles-methods
##' @aliases alleles-methods
##' @keywords methods
setGeneric("alleles<-",
           def=function(x, ..., value) {
             x@alleles <- value
             return(x)
           })


##' set phenotype variable in an annotSnpStats object
##'
##' @title phenotype
##' @param object object of class annotSnpStats
##' @param ... variable name that represents phenotype
##' @return no return value
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname phenotype-methods
##' @aliases phenotype-methods
##' @keywords methods
setGeneric("phenotype<-",
           def=function(x, ..., value) {
             x@phenotype <- value
             return(x)
           })

## setGeneric("align.alleles",
##            def=function(x, ref, ...) {
##              align.alleles(x, ref) } )

setGeneric("dimnames<-")

## setGeneric("ld",
##            def=function(x, ...) {
             
##            })
