##' Methods for aSnpMatrix and aXSnpMatrix objects
##'
##' These methods extend those that already exist for \code{matrix}
##' and \code{SnpMatrix} objects, so that the aligned snp and sample
##' summary data are correctly processed.
##'
##' @section Methods:
##' \describe{
##'
##' \item{\code{cbind2}}{binds two \code{aSnpMatrix} objects with the same samples, different SNPs}
##'
##' \item{\code{rbind2}}{binds two \code{aSnpMatrix} objects with the same SNPs, different samples}
##'
##' \item{\code{[}}{does the usual subsetting/reordering operations, but ties the SNP and sample support \code{data.frame}s in at the same time}
##'
##' }
##' 
##' @export
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases show,aSnpMatrix-method
##' @title Methods for aSnpMatrix objects
##' @param x aSnpMatrix or aXSnpMatrix object
##' @param object aSnpMatrix or aXSnpMatrix object
##' @param i row (sample) index 
##' @param j column (SNP) index 
##' @return object of same class as x
##' @author Chris Wallace
##' @examples
##' X <- new("aSnpMatrix")
##' X
##' ## load some example data from snpStats
##' data(for.exercise)
##' X <- new("aSnpMatrix",
##'          .Data = snps.10[1:10,1:5],
##'          snps=snp.support[1:5,],
##'          samples=subject.support[1:10,])
##' X
##' phenotype(X) <- "cc"
##' alleles(X) <- c("A1","A2")
##' X
##' summary(X)
setMethod("show",
          signature=c("aSnpMatrix"),
          function(object) {
            cat("Annotated SnpMatrix object with",nrow(object),"rows (samples) and",ncol(object),"columns (SNPs).\n")
            cat("phenotype:",if(length(object@phenotype)) {object@phenotype} else {"not set"},"\n")
            a <- suppressWarnings(alleles(object))
            cat("alleles:",if(length(a)) {paste(a,collapse="/")} else {"not set"},"\n")
          })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases [,aSnpMatrix,ANY,missing,missing-method
setMethod("[",
          signature=c(x="aSnpMatrix", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            if(is.logical(i))
              i <- which(i)
            new("aSnpMatrix",
                .Data=x@.Data[i,,drop=FALSE],
                snps=x@snps,
                samples=x@samples[i,,drop=FALSE],
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases [,aXSnpMatrix,ANY,missing,missing-method
##' @export
setMethod("[",
          signature=c(x="aXSnpMatrix", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            if(is.logical(i))
              i <- which(i)
            new("aXSnpMatrix",
                .Data=x@.Data[i,,drop=FALSE],
                snps=x@snps,
                samples=x@samples[i,,drop=FALSE],
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid[i])})
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases [,aSnpMatrix,missing,ANY,missing-method
setMethod("[",
          signature=c(x="aSnpMatrix", i="missing", j="ANY", drop="missing"),
          function(x, i, j) {
            if(is.logical(j))
              j <- which(j)
            new("aSnpMatrix",
##                .Data=new("SnpMatrix", matrix(as.raw(x@.Data),nrow=nrow(x@.Data), ncol=ncol(x@.Data))[,j,drop=FALSE] ),
                .Data=x@.Data[,j,drop=FALSE],
                snps=x@snps[j,,drop=FALSE],
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases [,aXSnpMatrix,missing,ANY,missing-method
setMethod("[",
          signature=c(x="aXSnpMatrix", i="missing", j="ANY", drop="missing"),
          function(x, i, j) {
            if(is.logical(j))
              j <- which(j)
            new("aXSnpMatrix",
##                .Data=new("SnpMatrix", matrix(as.raw(x@.Data),nrow=nrow(x@.Data), ncol=ncol(x@.Data))[,j,drop=FALSE] ),
                .Data=x@.Data[,j,drop=FALSE],
                snps=x@snps[j,,drop=FALSE],
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid)})

## system.time( kk<-new("SnpMatrix",
##                      matrix(as.raw(x@.Data),nrow=nrow(x@.Data), ncol=ncol(x@.Data),dimnames=dimnames(x@.Data))[,j,drop=FALSE]) )
## system.time( kk<-x@.Data[,j,drop=FALSE] )


##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases [,aSnpMatrix,ANY,ANY,missing-method
##' @aliases [,aSnpMatrix,missing,ANY,missing-method
setMethod("[",
          signature=c(x="aSnpMatrix", i="ANY", j="ANY", drop="missing"),
          function(x, i, j) {
            if(is.logical(i))
              i <- which(i)
     new("aSnpMatrix",
                .Data=x@.Data[i,j,drop=FALSE],
                snps=x@snps[j,,drop=FALSE],
                samples=x@samples[i,,drop=FALSE],
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases [,aXSnpMatrix,ANY,ANY,missing-method
##' @aliases [,aXSnpMatrix,missing,ANY,missing-method
setMethod("[",
          signature=c(x="aXSnpMatrix", i="ANY", j="ANY", drop="missing"),
          function(x, i, j) {
                if(is.logical(i))
              i <- which(i)
        new("aXSnpMatrix",
                .Data=x@.Data[i,j,drop=FALSE],
                snps=x@snps[j,],
                samples=x@samples[i,,drop=FALSE],
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid[i])})

fix.factors <- function(x) {
  if(!is.data.frame(x)) {
    warning("Object passed to fix,factors is not a data.frame. Returning unchanged.")
    return(x)
  }
  index.factors <- which(sapply(x,is.factor))
  if(length(index.factors))
    for(i in index.factors)
      x[,i] <- as.character(x[,i])
  return(x)
}
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' NB - takes allele labels from x, and assumes y matches.  NO CHECKS ARE MADE!  CAVEAT EMPTOR.
##' @aliases rbind2,aSnpMatrix,aSnpMatrix-method
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
                samples=rbind(fix.factors(x@samples[,samples.colmatch,drop=FALSE]),
                  fix.factors(y@samples[,samples.colmatch,drop=FALSE])),
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases rbind2,aXSnpMatrix,aXSnpMatrix-method
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

##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases cbind2,aSnpMatrix,aSnpMatrix-method
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
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases cbind2,aXSnpMatrix,aXSnpMatrix-method
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
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
setReplaceMethod("rownames",
                 signature=c(x="aSnpMatrix"),
                 function(x, value) {
                   if(length(value)!=nrow(x))
                     stop("must specify rownames with the same length as nrow(x)")
                   rownames(x@.Data) <- value
                   rownames(x@samples) <- value
                   new("aSnpMatrix",
                       .Data=x@.Data,
                       snps=x@snps,
                       samples=x@samples,
                       phenotype=x@phenotype,
                       alleles=x@alleles)
                 })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
setReplaceMethod("rownames",
                 signature=c(x="aXSnpMatrix"),
                 function(x, value) {
                   if(length(value)!=nrow(x))
                     stop("must specify rownames with the same length as nrow(x)")
                   rownames(x@.Data) <- value
                   rownames(x@samples) <- value
                   new("aXSnpMatrix",
                       .Data=x@.Data,
                       snps=x@snps,
                       samples=x@samples,
                       phenotype=x@phenotype,
                       alleles=x@alleles,
                       diploid=x@diploid)
                 })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
setReplaceMethod("colnames",
                 signature=c(x="aSnpMatrix"),
                 function(x, value) {
                   if(length(value)!=ncol(x))
                     stop("must specify colnames with the same length as ncol(x)")
                   colnames(x@.Data) <- value
                   rownames(x@snps) <- value
                   new("aSnpMatrix",
                       .Data=x@.Data,
                       snps=x@snps,
                       samples=x@samples,
                       phenotype=x@phenotype,
                       alleles=x@alleles)
                 })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
setReplaceMethod("colnames",
                 signature=c(x="aXSnpMatrix"),
                 function(x, value) {
                   if(length(value)!=ncol(x))
                     stop("must specify colnames with the same length as ncol(x)")
                   colnames(x@.Data) <- value
                   rownames(x@snps) <- value
                   new("aXSnpMatrix",
                       .Data=x@.Data,
                       snps=x@snps,
                       samples=x@samples,
                       phenotype=x@phenotype,
                       alleles=x@alleles,
                       diploid=x@diploid)
                 })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases dimnames<-,aSnpMatrix,ANY-method
setReplaceMethod("dimnames",
                 signature=c(x="aSnpMatrix"),
                 function(x, value) {
                   rownames(x) <- value[[1]]
                   colnames(x) <- value[[2]]
                   return(x) })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases dimnames<-,aXSnpMatrix,ANY-method
setReplaceMethod("dimnames",
                 signature=c(x="aXSnpMatrix"),
                 function(x, value) {
                   rownames(x) <- value[[1]]
                   colnames(x) <- value[[2]]
                   return(x) })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @param snps vector indexing SNPs whose alleles should be switched
##' @aliases switch.alleles,aSnpMatrix,ANY-method
setMethod("switch.alleles",
          signature=c(x="aSnpMatrix",snps="ANY"),
          function(x, snps) {
            x@.Data=switch.alleles(new("SnpMatrix",x@.Data), snps)
            anames <- alleles(x)
            x@snps[snps,anames] <- x@snps[snps,rev(anames)]
            return(x)
          })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
##' @aliases switch.alleles,aXSnpMatrix,ANY-method
setMethod("switch.alleles",
          signature=c(x="aXSnpMatrix",snps="ANY"),
          function(x, snps) {
            x@.Data=switch.alleles(new("SnpMatrix",x@.Data), snps)
            anames <- alleles(x)
            x@snps[snps,anames] <- x@snps[snps,rev(anames)]
            return(x)
          })


## setMethod("align.alleles",
##           signature=c(x="aSnpMatrix",ref="aSnpMatrix"),
##           align.alleles(x,ref))
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
setReplaceMethod("alleles",
          signature=c(x="aSnpMatrix"),
          function(x, value) {
            if(!(all(value %in% colnames(snps(x)))) || length(value)!=2)
              stop("alleles must be a vector of length 2 naming columns found in snps(x)")
            x@alleles <- value
            return(x)
          })
##' @rdname aSnpMatrix-methods
##' @name aSnpMatrixMethods
setReplaceMethod("phenotype",
          signature=c(x="aSnpMatrix"),
          function(x, value) {
            if(!(all(value %in% colnames(samples(x)))) || length(value)!=1)
              stop("phenotype must be the name of a column found in samples(x)")
           x@phenotype <- value
            return(x)
          })

##' @rdname add-methods
##' @name addMethods
##' @aliases add.snps,aSnpMatrix,aSnpMatrix-method
setMethod("add.snps",
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
                     function(x,y) {
                       tmp <- new("aSnpMatrix",
                                  .Data=matrix(as.raw("00"),nrow(x),ncol(y),dimnames=list(rownames(x),colnames(x))),
                                  snps=snps(y),
                                  samples=samples(x),
                                  phenotype=phenotype(x),
                                  alleles=alleles(y))
                       cbind2(x,tmp)
                     })

##' @rdname add-methods
##' @name addMethods
##' @aliases add.samples,aSnpMatrix,aSnpMatrix-method
setMethod("add.samples",
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
                     function(x,y) {
                       tmp <- new("aSnpMatrix",
                                  .Data=matrix(as.raw("00"),nrow(y),ncol(x),dimnames=list(rownames(y),colnames(x))),
                                  snps=snps(x),
                                  samples=samples(y),
                                  phenotype=phenotype(y),
                                  alleles=alleles(x))
                       rbind2(x,tmp)
                     })
