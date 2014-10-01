##' Binding
##' 
##' Bind two aSnpMatrix or aXSnpMatrix objects, or add empty
##' rows/columns to an existing object.
##'
##' \code{rbind2} and \code{cbind2} work as they do with SnpMatrix
##' objects
##' 
##' @param x object of class aSnpMatrix
##' @param y object of same class as \code{x}
##' @return object of same class as \code{x}
##' @author Chris Wallace
##' @export
##' @rdname binding
##' @keywords methods
##' @examples
##' ## load some example data from snpStats
##' data(for.exercise)
##' x <- new("aSnpMatrix",
##'          .Data = snps.10[1:10,1:5],
##'          snps=snp.support[1:5,],
##'          samples=subject.support[1:10,])
##' y <- new("aSnpMatrix",
##'          .Data = snps.10[1:10,11:15],
##'          snps=snp.support[11:15,],
##'          samples=subject.support[1:10,])
##' # bind x and y columnwise
##' z <- cbind2(x,y)
##' # add empty entries to x corresponding to the samples in x
##' # and the snps in y
##' z <- add.snps(x,y)
setMethod("rbind2",  ## bind samples
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
            if(!identical(colnames(x@.Data),colnames(y@.Data)))
               stop("SNP names mismatch")
            samples.colmatch <- intersect(colnames(x@samples),colnames(y@samples))
            ## missing alleles
            if(any(is.na(x@snps[,alleles(x)]))) {
              x.missing <- apply(is.na(x@snps[,alleles(x)]),1,any)
              y.missing <- apply(is.na(y@snps[,alleles(y)]),1,any)
              message("Missing alleles found.  Missingness table:")
              print(table(x.missing,y.missing))
              wh1 <- which(is.na(x@snps[,alleles(x)[1]]))
              wh2 <- which(is.na(x@snps[,alleles(x)[2]]))
              if(length(wh1))
                x@snps[wh1,alleles(x)[1]] <- y@snps[wh1,alleles(y)[1]] 
              if(length(wh2))
                x@snps[wh2,alleles(x)[2]] <- y@snps[wh2,alleles(y)[2]] 
              message("Missingness updated:")
              x.missing.update <- apply(is.na(x@snps[,alleles(x)]),1,any)
              print(table(x.missing,x.missing.update))
            }
            ## overlapping sample ids
            m <- match(rownames(x),rownames(y))
            if(any(!is.na(m))) {
              warning(sum(!is.na(m))," overlapping samples found - uniquifying sample names")
              new.ids <- make.unique(c(rownames(x)[!is.na(m)],rownames(y)[m[!is.na(m)]]))
              rownames(y)[m[!is.na(m)]] <- new.ids[-c(1:sum(!is.na(m)))]
            }
            new("aSnpMatrix",
                .Data=rbind2(as(x,"SnpMatrix"),
                  as(y,"SnpMatrix")),
                snps=x@snps,
                samples=rbind(x@samples[,samples.colmatch,drop=FALSE],
                  y@samples[,samples.colmatch,drop=FALSE]),
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname binding
##' @export
setMethod("rbind2", ## bind samples
          signature=c(x="aXSnpMatrix",y="aXSnpMatrix"),
          function(x,y) {
            if(!identical(colnames(x@.Data),colnames(y@.Data)))
              stop("SNP names mismatch")
            samples.colmatch <- intersect(colnames(x@samples),colnames(y@samples))
            new("aXSnpMatrix",
                .Data=rbind2(as(x,"XSnpMatrix"),
                  as(y,"XSnpMatrix")),
                snps=x@snps,
                samples=rbind(x@samples[,samples.colmatch,drop=FALSE],y@samples[,samples.colmatch,drop=FALSE]),
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=c(x@diploid,y@diploid))
          })

##' @rdname binding
##' @export
setMethod("cbind2", ## bind SNPs
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
            if(!identical(rownames(x@.Data),rownames(y@.Data)))
              stop("sample names mismatch")
            snps.colmatch <- intersect(colnames(x@snps),colnames(y@snps))
            new("aSnpMatrix",
                .Data=cbind2(as(x,"SnpMatrix"),
                  as(y,"SnpMatrix")),
                snps=rbind(x@snps[,snps.colmatch,drop=FALSE],y@snps[,snps.colmatch,drop=FALSE]),
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname binding
##' @export
setMethod("cbind2", ## bind SNPs
          signature=c(x="aXSnpMatrix",y="aXSnpMatrix"),
          function(x,y) {
            if(!identical(rownames(x@.Data),rownames(y@.Data)))
              stop("sample names mismatch")
            if(!identical(x@diploid, y@diploid))
              stop("sample diploid status mismatch")
            snps.colmatch <- intersect(colnames(x@snps),colnames(y@snps))
            new("aXSnpMatrix",
                .Data=cbind2(as(x,"XSnpMatrix"),
                  as(y,"XSnpMatrix")),
                snps=rbind(x@snps[,snps.colmatch,drop=FALSE],y@snps[,snps.colmatch,drop=FALSE]),
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid)})


##' @export
##' @rdname binding
##' @details \code{add.snps()}  adds empty entries to x for SNPs found in y
setGeneric("add.snps",
           signature=c("x","y"),
           def=function(x,y) standardGeneric("add.snps"))
##' @export
##' @rdname binding
setMethod("add.snps",
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
               tmp <- new("aSnpMatrix",
                          .Data=matrix(as.raw("00"),nrow(x),ncol(y),
                              dimnames=list(rownames(x),colnames(y))),
                          snps=snps(y),
                          samples=samples(x),
                          phenotype=phenotype(x),
                          alleles=alleles(y))
               cbind2(x,tmp)
           })

##' @export
##' @rdname binding
##' @details \code{add.samples()}  adds empty entries to x for Samples found in y
setGeneric("add.samples",
           signature=c("x","y"),
           def=function(x,y) standardGeneric("add.samples"))
##' @rdname binding
##' @export
setMethod("add.samples",signature=c(x="aSnpMatrix",y="aSnpMatrix"),
           function(x,y) {
                       tmp <- new("aSnpMatrix",
                                  .Data=matrix(as.raw("00"),nrow(y),ncol(x),dimnames=list(rownames(y),colnames(x))),
                                  snps=snps(x),
                                  samples=samples(y),
                                  phenotype=phenotype(y),
                                  alleles=alleles(x))
                       rbind2(x,tmp)
                     })