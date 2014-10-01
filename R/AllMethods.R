##' Methods for aSnpMatrix and aXSnpMatrix objects
##'
##' These methods extend those that already exist for \code{matrix}
##' and \code{SnpMatrix} objects, so that the aligned snp and sample
##' summary data are correctly processed.
##'
##' \code{cbind2} binds two \code{aSnpMatrix} objects with the same samples, different SNPs
##'
##' \code{rbind2} binds two \code{aSnpMatrix} objects with the same SNPs, different samples
##'
##' \code{[} does the usual subsetting/reordering operations, but ties the SNP and sample support \code{data.frame}s in at the same time
##'
##' @export
##' @rdname methods
##' @param object aSnpMatrix or aXSnpMatrix object
##' @param x aSnpMatrix or aXSnpMatrix object
##' @param i row (sample) index 
##' @param j column (SNP) index 
##' @param value replacement value
##' @param snps index vector of SNPs for which alleles should be switched
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
setMethod("show",
          signature=c("aSnpMatrix"),
          function(object) {
            cat("Annotated SnpMatrix object with",nrow(object),"rows (samples) and",ncol(object),"columns (SNPs).\n")
            cat("phenotype:",if(length(object@phenotype)) {object@phenotype} else {"not set"},"\n")
            a <- suppressWarnings(alleles(object))
            cat("alleles:",if(length(a)) {paste(a,collapse="/")} else {"not set"},"\n")
          })

##' @rdname methods
##' @export
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
##' @rdname methods
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

##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
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

##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
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
##' @rdname methods
##' @export
setReplaceMethod("dimnames",
                 signature=c(x="aSnpMatrix"),
                 function(x, value) {
                   rownames(x) <- value[[1]]
                   colnames(x) <- value[[2]]
                   return(x) })
##' @rdname methods
##' @export
setReplaceMethod("dimnames",
                 signature=c(x="aXSnpMatrix"),
                 function(x, value) {
                   rownames(x) <- value[[1]]
                   colnames(x) <- value[[2]]
                   return(x) })
##' @rdname methods
##' @export
setMethod("switch.alleles",
          signature=c(x="aSnpMatrix",snps="ANY"),
          function(x, snps) {
            x@.Data=switch.alleles(new("SnpMatrix",x@.Data), snps)
            anames <- alleles(x)
            x@snps[snps,anames] <- x@snps[snps,rev(anames)]
            return(x)
          })
##' @rdname methods
##' @export
setMethod("switch.alleles",
          signature=c(x="aXSnpMatrix",snps="ANY"),
          function(x, snps) {
            x@.Data=switch.alleles(new("SnpMatrix",x@.Data), snps)
            anames <- alleles(x)
            x@snps[snps,anames] <- x@snps[snps,rev(anames)]
            return(x)
          })

