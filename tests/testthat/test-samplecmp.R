data(for.exercise, package="snpStats")
n <- 500
m <- 20000

context("dups")

x <- snps.10[1:n,1:m]
y <- snps.10[1:(2*n),1:m]
rownames(x) <- paste("x",1:n,sep="")
rownames(y) <- paste("y",1:(2*n),sep="")

test_that("dups works", {
  kk <- dups(x,y,tol=10)
  expect_true(nrow(kk) == n)
})
          
## system.time(kk<-dups(x,y))
## system.time(ll<-dups.c(x,y))

## hist(kk,breaks=m)
## hist(diag(kk),breaks=m,xlim=c(0,m))
## hist(kk[upper.tri(kk)],breaks=m,xlim=c(0,m))
## hist(kk[lower.tri(kk)],breaks=m,xlim=c(0,m))
