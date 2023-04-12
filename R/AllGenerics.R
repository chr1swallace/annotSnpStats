
##' Accessors
##' 
##' get/set variables in an annotSnpStats object
##'
##' \code{alleles()} gets the two element character vector which
##' names the columns in the \code{snps} slot which records the
##' ordering of base alleles for that SNP.
##'
##' @param x object of class annotSnpStats
##' @param value replacement value
##' @return the column names requested
##' @author Chris Wallace
##' @export
##' @rdname accessors
##' @keywords methods
##' @examples
##' X <- example.data(10,5)
##' phenotype(X)
##' phenotype(X) <- "cc"
##' alleles(X)
##' snps(X)
##' samples(X)
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
##' @rdname accessors
##' @details
##' For \code{alleles<-}, value is a character vector of
##' length 2, giving the columns names of the data.frame in the
##' \code{snps} slot corresponding to alelles 1 and 2.
setGeneric("alleles<-",
           def=function(x, value) {
               if(!(all(value %in% colnames(snps(x)))) || length(value)!=2)
                   stop("alleles must be a vector of length 2 naming columns found in snps(x)")
               x@alleles <- value
               return(x)
           })


##' @export
##' @rdname accessors
##' @details
##' \code{phenotype} gets/sets the column name of the data.frame in the
##' \code{samples} slots which encodes the phenotype.
setGeneric("phenotype",
           def=function(x) {
                     x@phenotype
           })
##' @export
##' @rdname accessors
##' @aliases phenotype<-,aSnpMatrix
##' @aliases phenotype<-,aXSnpMatrix
setGeneric("phenotype<-",
           def=function(x, value) {
             if(!(value %in% colnames(samples(x))) || length(value)!=1)
               stop("phenotype must be name of a column in samples(x)")
             x@phenotype <- value
             return(x)
           })


##' @export
##' @rdname accessors
##' @details \code{snps()} extracts the snps \code{data.frame} from an
##' object of class aSnpMatrix
setGeneric("snps",
           def=function(x) {
             return(x@snps)
           })
##' @export
##' @rdname accessors
setGeneric("snps<-",
           def=function(x,value) {
               if(!is.data.frame(value))
                   stop("snps must be a data.frame")
               x@snps <- value
               return(object)
           })
##' @export
##' @rdname accessors
##' @details \code{samples()} extracts the samples \code{data.frame}
##' from an object of class aSnpMatrix
setGeneric("samples",
           def=function(x) {
             return(x@samples)
           })

##' @export
##' @rdname accessors
##' @details \code{samples()} extracts the samples \code{data.frame}
##' from an object of class aSnpMatrix
setGeneric("samples<-",
           def=function(x,value) {
             if(!is.data.frame(value))
               stop("samples must be a data.frame")
             if(!identical(rownames(x),rownames(value)))
               stop("samples must have identical nrow and rownames to x")
             x@samples <- value
             return(x)
           })




## ##' set rownames and colnames simultaneously
## ##'
## ##' @export
## ##' @param x aSnpMatrix or aXSnpMatrix object
## ##' @param value list of two character vectors giving row and sample names
## setGeneric("dimnames<-")
