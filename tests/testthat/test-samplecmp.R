data(for.exercise, package="snpStats")
n <- 700
m <- 20000

context("dups")

x <- snps.10[c(1:n,1:n),1:m]
y <- snps.10[c(1:1000,1:1000),1:m]
rownames(x) <- paste("x",make.unique(rownames(x)),sep="")
rownames(y) <- paste("y",make.unique(rownames(y)),sep="")

test_that("dups works", {
  system.time(kk.slow <- dups(x,y,tol=10,stopatone=FALSE))
  system.time(kk.quick <- dups(x,y,tol=10,stopatone=TRUE))
  expect_true(nrow(kk) == n)
})
          
## system.time(kk<-dups(x,y))
## system.time(ll<-dups.c(x,y))

## hist(kk,breaks=m)
## hist(diag(kk),breaks=m,xlim=c(0,m))
## hist(kk[upper.tri(kk)],breaks=m,xlim=c(0,m))
## hist(kk[lower.tri(kk)],breaks=m,xlim=c(0,m))
