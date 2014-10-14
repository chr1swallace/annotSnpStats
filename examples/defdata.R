data(for.exercise, package="snpStats")
n <- 10
m <- 100
rownames(Asnps) <- colnames(Autosomes)
Asnps$A1 <- "A"
Asnps$A2 <- "G"
x <- new("aSnpMatrix",
            .Data=Autosomes[1:m,1:n],
            snps=Asnps[1:n,,drop=FALSE],
            samples=subject.data[1:m,],
            phenotype="cc")
y <- new("aSnpMatrix",
                 .Data=Autosomes[(m+1):(2*m),1:n],
                 snps=Asnps[1:n,,drop=FALSE],
                 samples=subject.data[(m+1):(2*m),],
            phenotype="cc")
