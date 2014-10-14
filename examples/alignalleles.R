# switch columns 6:10
y.switched <- switch.alleles(y,6:10)
# automatically switch back
y.aligned <- align.alleles(y.switched,x,do.plot=TRUE)
## check by comparing counted allele frequencies for SNPs 1:5 and 6:10
cbind(x=col.summary(x)[,"RAF"],
      y=col.summary(y)[,"RAF"], 
      y.switched=col.summary(y.switched)[,"RAF"],
      y.aligned=col.summary(y.aligned)[,"RAF"])
     
