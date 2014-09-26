library(snpStats)
data(testdata, package="snpStats")
n <- 10
m <- 400

context("impute")
rownames(Asnps) <- colnames(Autosomes)
rownames(Xsnps) <- colnames(Xchromosome)

i1 <- 1:400
j1 <- 101:120
toimp <- sample
x <- Autosomes[i1,j1]
x2<-x
x2[cbind(1:10,1:10)] <- as.raw("00")
x2 <- new("aSnpMatrix",
          .Data=x2,
          snps=Asnps[j1,,drop=FALSE],
          samples=subject.data)

y <- impute.missing(x2) 
as(y,"numeric")[cbind(1:10,1:10)]
as(x,"numeric")[cbind(1:10,1:10)]
