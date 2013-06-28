data(for.exercise, package="snpStats")
n <- 20
m <- 100

s10 <- snps.10[1:m,1:n]

## initialize

as10 <- new("aSnpMatrix",
                 .Data=s10,
                 snps=snp.support[1:n,],
                 samples=subject.support[1:m,])
as10.2 <- new("aSnpMatrix",
                 .Data=snps.10[1:m, n+(1:n)],
                 snps=snp.support[n + (1:n),],
                 samples=subject.support[1:m,])
as10.3 <- new("aSnpMatrix",
                 .Data=snps.10[m+(1:m), 1:n],
                 snps=snp.support[1:n,],
                 samples=subject.support[m+(1:m),])
test_that("initialize works", {
  expect_that(as10, is_a("aSnpMatrix"))
  expect_that(new("aSnpMatrix",
                  .Data=s10,
                  snps=snp.support,
                  samples=subject.support), throws_error())
  expect_that(new("aSnpMatrix",
                  .Data=s10,
                  snps=snp.support[1:n,],
                  samples=subject.support[1:n,]), throws_error())
})

## subsetting
test_that("subsetting works", {
  expect_that(nrow(as10[1:10,]), equals(10))
  expect_that(ncol(as10[,1:10]), equals(10))
})

## binding
test_that("binding works", {
  expect_that(ncol(cbind2(as10,as10.2)), equals(n * 2))
  expect_that(nrow(rbind2(as10,as10.3)), equals(m * 2))
  expect_that(cbind2(as10,as10.3), throws_error())
  expect_that(rbind2(as10,as10.2), throws_error())
})

          
          
