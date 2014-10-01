data(testdata, package="snpStats")
n <- 20
m <- 100

i1 <- 1:m
i2 <- m+(1:m)
j1 <- 1:n
j2 <- n+(1:n)


## initialize
rownames(Asnps) <- colnames(Autosomes)
rownames(Xsnps) <- colnames(Xchromosome)
Asnps$a1 <- "A"
Asnps$a2 <- "B"
Xsnps$a1 <- "A"
Xsnps$a2 <- "B"


as11 <- new("aSnpMatrix",
            .Data=Autosomes[i1,j1],
            snps=Asnps[j1,,drop=FALSE],
            samples=subject.data[i1,],
            phenotype="cc")
as12 <- new("aSnpMatrix",
                 .Data=Autosomes[i1,j2],
                 snps=Asnps[j2,,drop=FALSE],
                 samples=subject.data[i1,],
            phenotype="cc")
as21 <- new("aSnpMatrix",
                 .Data=Autosomes[i2,j1],
                 snps=Asnps[j1,,drop=FALSE],
                 samples=subject.data[i2,],
            phenotype="cc")
xs11 <- new("aXSnpMatrix",
                 .Data=Xchromosome[i1,j1],
                 snps=Xsnps[j1,,drop=FALSE],
                 samples=subject.data[i1,],
            phenotype="cc")
xs12 <- new("aXSnpMatrix",
                 .Data=Xchromosome[i1,j2],
                 snps=Xsnps[j2,,drop=FALSE],
                 samples=subject.data[i1,],
            phenotype="cc")
xs21 <- new("aXSnpMatrix",
                 .Data=Xchromosome[i2,j1],
                 snps=Xsnps[j1,,drop=FALSE],
                 samples=subject.data[i2,],
            phenotype="cc")


context("initialize")
test_that("initialize works", {
  expect_that(as11, is_a("aSnpMatrix"))
  expect_that(xs11, is_a("aXSnpMatrix"))
  expect_that(new("aSnpMatrix",
                  .Data=as11@.Data,
                  snps=as12@snps,
                  samples=as11@samples), throws_error())
  expect_that(new("aXSnpMatrix",
                  .Data=xs11@.Data,
                  snps=xs12@snps,
                  samples=xs11@samples), throws_error())
  expect_that(new("aSnpMatrix",
                  .Data=as11@.Data,
                  snps=as11@snps,
                  samples=as21@samples), throws_error())
  expect_that(new("aXSnpMatrix",
                  .Data=xs11@.Data,
                  snps=xs11@snps,
                  samples=xs21@samples), throws_error())
})

## subsetting
context("subsetting")
test_that("subsetting works", {
  expect_that(as11[1:10,], is_a("aSnpMatrix"))
  expect_that(xs11[1:10,], is_a("aXSnpMatrix"))
  expect_that(nrow(as11[1:10,]), equals(10))
  expect_that(ncol(as11[,1:10]), equals(10))
  expect_that(nrow(xs11[1:10,]), equals(10))
  expect_that(ncol(xs11[,1:10]), equals(10))
})

## binding
context("binding")
test_that("binding works", {
  ac <- cbind2(as11,as12)
  ar <- rbind2(as11,as21)
  xc <- cbind2(xs11,xs12)
  xr <- rbind2(xs11,xs21)

  expect_that(ac, is_a("aSnpMatrix"))
  expect_that(ar, is_a("aSnpMatrix")) 
  expect_that(xc, is_a("aXSnpMatrix"))
  expect_that(xr, is_a("aXSnpMatrix"))
  expect_that(ncol(ac), equals(n * 2))
  expect_that(nrow(ar), equals(m * 2))
  expect_that(ncol(xc), equals(n * 2))
  expect_that(nrow(xr), equals(m * 2))
 
  expect_that(cbind2(as11,as21), throws_error())
  expect_that(rbind2(as11,as12), throws_error())
  expect_that(cbind2(xs11,xs21), throws_error())
  expect_that(rbind2(xs11,xs12), throws_error())
  expect_that(cbind2(as11,as11), throws_error())
  expect_that(cbind2(xs11,xs11), throws_error())
  expect_warning(rbind2(as11,as11))
  expect_that(rbind2(xs11,xs11), throws_error())
})

## rownames
context("dimnames")
test_that("dimnames/rownames/colnames", {
  nm <- paste("s",1:m,sep="")
  nn <- paste("s",1:n,sep="")
  for(tmp in list(as11, xs11)) {
    dimnames(tmp) <- list(nm,nn)
    expect_identical(rownames(tmp@.Data), nm)
    expect_identical(rownames(tmp@samples), nm)
    expect_identical(colnames(tmp@.Data), nn)
    expect_identical(rownames(tmp@snps), nn)
  }
})

## conversion
context("conversion")
test_that("conversion", {
  expect_that(as(as11,"SnpMatrix"), is_a("SnpMatrix"))
  expect_that(as(xs11,"XSnpMatrix"), is_a("XSnpMatrix"))
  expect_warning(as(as11,"XSnpMatrix"))
  expect_that(as(xs11,"SnpMatrix"), is_a("SnpMatrix"))
})
          
