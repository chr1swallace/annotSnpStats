##' Methods for aSnpMatrix objects
##'
##' Most methods simply extend methods that already exist for \code{SnpMatrix} objects.
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
##' @docType methods
##' @rdname aSnpMatrix-methods
##' @aliases [,aSnpMatrix,ANY,missing,missing-method
##' @title Methods for aSnpMatrix objects
##' @param x aSnpMatrix object
##' @param y aSnpMatrix object
##' @param i row (sample) index 
##' @param j column (SNP) index 
##' @return aSnpMatrix object
##' @author Chris Wallace
setMethod("[",
          signature=c(x="aSnpMatrix", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            new("aSnpMatrix",
                .Data=x@.Data[i,],
                snps=x@snps,
                samples=x@samples[i,],
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @aliases [,aXSnpMatrix,missing,ANY,missing-method
setMethod("[",
          signature=c(x="aXSnpMatrix", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            new("aXSnpMatrix",
                .Data=x@.Data[i,],
                snps=x@snps,
                samples=x@samples[i,],
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid[i])})
##' @name [
##' @rdname aSnpMatrix-methods
##' @aliases [,aSnpMatrix,missing,ANY,missing-method
##' @docType methods
setMethod("[",
          signature=c(x="aSnpMatrix", i="missing", j="ANY", drop="missing"),
          function(x, i, j) {
            new("aSnpMatrix",
##                .Data=new("SnpMatrix", matrix(as.raw(x@.Data),nrow=nrow(x@.Data), ncol=ncol(x@.Data))[,j,drop=FALSE] ),
                .Data=x@.Data[,j],
                snps=x@snps[j,,drop=FALSE],
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @name [
##' @rdname aSnpMatrix-methods
##' @aliases [,aXSnpMatrix,missing,ANY,missing-method
##' @docType methods
setMethod("[",
          signature=c(x="aXSnpMatrix", i="missing", j="ANY", drop="missing"),
          function(x, i, j) {
            new("aXSnpMatrix",
##                .Data=new("SnpMatrix", matrix(as.raw(x@.Data),nrow=nrow(x@.Data), ncol=ncol(x@.Data))[,j,drop=FALSE] ),
                .Data=x@.Data[,j],
                snps=x@snps[j,,drop=FALSE],
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid)})

## system.time( kk<-new("SnpMatrix",
##                      matrix(as.raw(x@.Data),nrow=nrow(x@.Data), ncol=ncol(x@.Data),dimnames=dimnames(x@.Data))[,j,drop=FALSE]) )
## system.time( kk<-x@.Data[,j,drop=FALSE] )


##' @rdname aSnpMatrix-methods
##' @aliases [,aSnpMatrix,ANY,ANY,missing-method
setMethod("[",
          signature=c(x="aSnpMatrix", i="ANY", j="ANY", drop="missing"),
          function(x, i, j) {
            new("aSnpMatrix",
                .Data=x@.Data[i,j],
                snps=x@snps[j,,drop=FALSE],
                samples=x@samples[i,,drop=FALSE],
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @aliases [,aXSnpMatrix,ANY,ANY,missing-method
setMethod("[",
          signature=c(x="aXSnpMatrix", i="ANY", j="ANY", drop="missing"),
          function(x, i, j) {
            new("aXSnpMatrix",
                .Data=x@.Data[i,j],
                snps=x@snps[j,],
                samples=x@samples[i,],
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid[i])})

##' @rdname aSnpMatrix-methods
##' @aliases rbind2,aSnpMatrix,aSnpMatrix
setMethod("rbind2",  ## bind samples
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
            if(!identical(colnames(x@.Data),colnames(y@.Data)))
               stop("SNP names mismatch")
            samples.colmatch <- intersect(colnames(x@samples),colnames(y@samples))
            new("aSnpMatrix",
                .Data=rbind2(x@.Data,y@.Data),
                snps=x@snps,
                samples=rbind(x@samples[,samples.colmatch,drop=FALSE],y@samples[,samples.colmatch,drop=FALSE]),
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aSnpMatrix-methods
##' @aliases rbind2,aXSnpMatrix,aXSnpMatrix
setMethod("rbind2", ## bind samples
          signature=c(x="aXSnpMatrix",y="aXSnpMatrix"),
          function(x,y) {
            if(!identical(colnames(x@.Data),colnames(y@.Data)))
              stop("SNP names mismatch")
             samples.colmatch <- intersect(colnames(x@samples),colnames(y@samples))
            new("aXSnpMatrix",
                .Data=rbind2(x@.Data,y@.Data),
                snps=x@snps,
                samples=rbind(x@samples[,samples.colmatch,drop=FALSE],y@samples[,samples.colmatch,drop=FALSE]),
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=c(x@diploid,y@diploid)) })

##' @rdname aSnpMatrix-methods
##' @aliases cbind2,aSnpMatrix,aSnpMatrix
setMethod("cbind2", ## bind SNPs
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
            if(!identical(rownames(x@.Data),rownames(y@.Data)))
              stop("sample names mismatch")
            snps.colmatch <- intersect(colnames(x@snps),colnames(y@snps))
            new("aSnpMatrix",
                .Data=cbind2(x@.Data,y@.Data),
                snps=rbind(x@snps[,snps.colmatch,drop=FALSE],y@snps[,snps.colmatch,drop=FALSE]),
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles)})
##' @rdname aXSnpMatrix-methods
##' @aliases cbind2,aXSnpMatrix,aXSnpMatrix
setMethod("cbind2", ## bind SNPs
          signature=c(x="aXSnpMatrix",y="aXSnpMatrix"),
          function(x,y) {
            if(!identical(rownames(x@.Data),rownames(y@.Data)))
              stop("sample names mismatch")
            if(!identical(x@diploid, y@diploid))
              stop("sample diploid status mismatch")
            snps.colmatch <- intersect(colnames(x@snps),colnames(y@snps))
            new("aXSnpMatrix",
                .Data=cbind2(x@.Data,y@.Data),
                snps=rbind(x@snps[,snps.colmatch,drop=FALSE],y@snps[,snps.colmatch,drop=FALSE]),
                samples=x@samples,
                phenotype=x@phenotype,
                alleles=x@alleles,
                diploid=x@diploid)})

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
setReplaceMethod("dimnames",
                 signature=c(x="aSnpMatrix"),
                 function(x, value) {
                   rownames(x) <- value[[1]]
                   colnames(x) <- value[[2]]
                   return(x) })
setReplaceMethod("dimnames",
                 signature=c(x="aXSnpMatrix"),
                 function(x, value) {
                   rownames(x) <- value[[1]]
                   colnames(x) <- value[[2]]
                   return(x) })
setMethod("switch.alleles",
          signature=c(x="aSnpMatrix",snps="ANY"),
          function(x, snps) {
            x@.Data=switch.alleles(new("SnpMatrix",x@.Data), snps)
            anames <- alleles(x)
            x@snps[snps,anames] <- x@snps[snps,rev(anames)]
            return(x)
          })
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
## setReplaceMethod("alleles",
##           signature=c(x="aSnpMatrix"),
##           function(x, value) {
##             x@alleles <- value
##             return(x)
##           })
## setReplaceMethod("phenotype",
##           signature=c(x="aSnpMatrix"),
##           function(x, value) {
##             x@phenotype <- value
##             return(x)
##           })
