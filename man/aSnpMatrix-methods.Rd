\docType{methods}
\name{[}
\alias{[}
\alias{[,aSnpMatrix,ANY,ANY,missing-method}
\alias{[,aSnpMatrix,ANY,missing,missing-method}
\alias{[,aSnpMatrix,missing,ANY,missing-method}
\alias{cbind2}
\alias{cbind2,aSnpMatrix,aSnpMatrix}
\alias{rbind2}
\alias{rbind2,aSnpMatrix,aSnpMatrix}
\title{Methods for aSnpMatrix objects}
\arguments{
  \item{x}{aSnpMatrix object}

  \item{y}{aSnpMatrix object}

  \item{i}{row (sample) index}

  \item{j}{column (SNP) index}
}
\value{
  aSnpMatrix object
}
\description{
  Methods for aSnpMatrix objects
}
\details{
  Most methods simply extend methods that already exist for
  \code{SnpMatrix} objects.
}
\section{Methods}{
  \describe{

  \item{\code{cbind2}}{binds two \code{aSnpMatrix} objects
  with the same samples, different SNPs}

  \item{\code{rbind2}}{binds two \code{aSnpMatrix} objects
  with the same SNPs, different samples}

  \item{\code{[}}{does the usual subsetting/reordering
  operations, but ties the SNP and sample support
  \code{data.frame}s in at the same time}

  }
}
\author{
  Chris Wallace
}

