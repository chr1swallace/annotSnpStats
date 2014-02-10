\docType{methods}
\name{aSnpMatrix-methods}
\alias{[}
\alias{[,aSnpMatrix,ANY,ANY,missing-method}
\alias{[,aSnpMatrix,ANY,missing,missing-method}
\alias{[,aSnpMatrix,missing,ANY,missing-method}
\alias{[,aXSnpMatrix,ANY,ANY,missing-method}
\alias{[,aXSnpMatrix,ANY,missing,missing-method}
\alias{[,aXSnpMatrix,missing,ANY,missing-method}
\alias{aSnpMatrix-methods}
\alias{cbind2,aSnpMatrix,aSnpMatrix-method}
\alias{cbind2,aXSnpMatrix,aXSnpMatrix-method}
\alias{colnames<-,aSnpMatrix-method}
\alias{colnames<-,aXSnpMatrix-method}
\alias{dimnames<-,aSnpMatrix-method}
\alias{dimnames<-,aXSnpMatrix-method}
\alias{rbind2,aSnpMatrix,aSnpMatrix-method}
\alias{rbind2,aXSnpMatrix,aXSnpMatrix-method}
\alias{rownames,aSnpMatrix-method}
\alias{rownames<-,aXSnpMatrix-method}
\alias{switch.alleles,aSnpMatrix,ANY}
\alias{switch.alleles,aXSnpMatrix,ANY}
\title{Methods for aSnpMatrix objects}
\usage{
\S4method{[}{aSnpMatrix,ANY,missing,missing}(x, i)

\S4method{[}{aXSnpMatrix,ANY,missing,missing}(x, i)

\S4method{[}{aSnpMatrix,missing,ANY,missing}(x, i, j)

\S4method{[}{aXSnpMatrix,missing,ANY,missing}(x, i, j)

\S4method{[}{aSnpMatrix,ANY,ANY,missing}(x, i, j)

\S4method{[}{aXSnpMatrix,ANY,ANY,missing}(x, i, j)

\S4method{rbind2}{aSnpMatrix,aSnpMatrix}(x, y)

\S4method{rbind2}{aXSnpMatrix,aXSnpMatrix}(x, y)

\S4method{cbind2}{aSnpMatrix,aSnpMatrix}(x, y)

\S4method{cbind2}{aXSnpMatrix,aXSnpMatrix}(x, y)

\S4method{switch.alleles}{aSnpMatrix,ANY}(x, snps)

\S4method{switch.alleles}{aXSnpMatrix,ANY}(x, snps)
}
\arguments{
  \item{x}{aSnpMatrix object}

  \item{y}{aSnpMatrix object}

  \item{i}{row (sample) index}

  \item{j}{column (SNP) index}

  \item{snps}{A vector of type integer, character or
  logical specifying the SNP to have its alleles switched}
}
\value{
aSnpMatrix object
}
\description{
Methods for aSnpMatrix objects
}
\details{
These methods extend those that already exist for
\code{matrix} and \code{SnpMatrix} objects, so that the
aligned snp and sample summary data are correctly
processed.
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

