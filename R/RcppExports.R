# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cdupsw <- function(X, maxDiff, Rtype, Rquick, null, het) {
    .Call('_annotSnpStats_cdupsw', PACKAGE = 'annotSnpStats', X, maxDiff, Rtype, Rquick, null, het)
}

cdups <- function(X, Y, maxDiff, Rtype, Rquick, null, het) {
    .Call('_annotSnpStats_cdups', PACKAGE = 'annotSnpStats', X, Y, maxDiff, Rtype, Rquick, null, het)
}

asw <- function(x, w) {
    .Call('_annotSnpStats_asw', PACKAGE = 'annotSnpStats', x, w)
}

