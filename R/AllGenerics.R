setGeneric("alleles<-",
           def=function(x, ..., value) {
             x@alleles <- value
             return(x)
           })
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

setGeneric("phenotype<-",
           def=function(x, ..., value) {
             x@phenotype <- value
             return(x)
           })

## setGeneric("align.alleles",
##            def=function(x, ref, ...) {
##              align.alleles(x, ref) } )

