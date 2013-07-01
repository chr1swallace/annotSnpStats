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
                samples=x@samples[i,])})
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
                snps=x@snps[j,],
                samples=x@samples)})

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
                snps=x@snps[j,],
                samples=x@samples[i,])})

##' @rdname aSnpMatrix-methods
##' @aliases rbind2,aSnpMatrix,aSnpMatrix
setMethod("rbind2",
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
            samples.colmatch <- intersect(colnames(x@samples),colnames(y@samples))
            new("aSnpMatrix",
                .Data=rbind2(x@.Data,y@.Data),
                snps=x@snps,
                samples=rbind(x@samples[,samples.colmatch],y@samples[,samples.colmatch]))})
##' @rdname aSnpMatrix-methods
##' @aliases cbind2,aSnpMatrix,aSnpMatrix
setMethod("cbind2",
          signature=c(x="aSnpMatrix",y="aSnpMatrix"),
          function(x,y) {
            snps.colmatch <- intersect(colnames(x@snps),colnames(y@snps))
            new("aSnpMatrix",
                .Data=cbind2(x@.Data,y@.Data),
                snps=rbind(x@snps[,snps.colmatch],y@snps[,snps.colmatch]),
                samples=x@samples)})
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
                       samples=x@samples)
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
                       samples=x@samples)
                 })

setMethod("switch.alleles",
          signature=c(x="aSnpMatrix",snps="ANY"),
          function(x, snps) {
            x@.Data=switch.alleles(new("SnpMatrix",x@.Data), snps)
            alleles <- x@snps[,c("allele.1","allele.2")]
            x@snps[sw,c("allele.1","allele.2")] <- x@snps[sw,c("allele.2","allele.1")]
            return(x)
          })
          
