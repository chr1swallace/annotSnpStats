\name{dups}
\alias{dups}
\title{Find duplicate samples}
\usage{
  dups(x, y, tol = ncol(x)/50, type = c("hethom", "all"))
}
\arguments{
  \item{x}{aSnpStats object}

  \item{y}{aSnpStats object}

  \item{tol}{maximum number of mismatched genotypes allowed
  for duplicate samples}

  \item{type}{by default, dups compares only homs vs hets,
  to allow for differently labelled alleles.  Set
  type="all" to allow the two kinds of homozygote genotypes
  to count as a mismatch}
}
\value{
  a matrix, with four columns: index of dup in x, index of
  dup in y, number of mismatches, number of comparisons
}
\description{
  Find indices of possible sample duplications between two
  aSnpStats objects
}
\details{
  Each pair of samples from x and y are compared in turn.
  If the number of mismatched and non-missing genotypes
  exceeds tol, the pair are assumed to be non-duplicates,
  and counting proceeds to the next pair.  If the total
  number of mismatched and non-missing genotypes is <tol,
  then the indices of the sample pair are stored, and
  returned together with the number of mismatches and the
  number of non-missing genotypes compared.
}
\author{
  Chris Wallace
}

