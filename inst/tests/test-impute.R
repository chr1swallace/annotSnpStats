data(for.exercise, package="snpStats")
n <- 500
m <- 20000

context("impute")

x <- snps.10[1:n,1:m]
y <- snps.10[1:(2*n),1:m]
