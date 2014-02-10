\name{snp.match}
\alias{snp.match}
\title{snp.match}
\usage{
snp.match(x, y, col.x = 0, col.y = 0)
}
\arguments{
  \item{x}{aSnpMatrix object}

  \item{y}{aSnpMatrix object}

  \item{col.x}{optional, column identifier for
  \code{x@snps} for matching}

  \item{col.y}{optional, column identifier for
  \code{y@snps} for matching}
}
\value{
a named list of length two with elements \describe{

\item{x}{indices of matching SNPs in x}

\item{y}{indices of matching SNPs in y}

}
}
\description{
Match snps in two aSnpMatrix objects
}
\details{
Extract indices of overlapping SNPs in two aSnpMatrix
objects.  By default, match is on \code{rownames(x@snps)},
\code{rownames(y@snps)}.
}
\author{
Chris Wallace
}

