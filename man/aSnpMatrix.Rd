\name{aSnpMatrix}
\alias{aSnpMatrix}
\alias{aSnpMatrix-class}
\title{Annotated SnpMatrix class}
\description{
  SnpMatrix objects, tied to sample and SNP support objects
}
\details{
  This is a small class to allow subsetting and binding
  operations to be applied to SnpMatrix objects together
  with their SNP and sample support objects in single
  operations.  It means the different objects are assured
  to always line up.
}
\section{Slots}{
  \describe{ \item{\code{.Data}:}{\code{SnpMatrix} object}
  \item{\code{snps}:}{\code{data.frame} with rownames
  matching colnames of \code{.Data}}
  \item{\code{samples}:}{\code{data.frame} with rownames
  matching rownames of \code{.Data}} }
}
\author{
  Chris Wallace
}

