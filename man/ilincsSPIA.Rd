% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ilincsSPIA.R
\name{ilincsSPIA}
\alias{ilincsSPIA}
\title{A multicore method to run SPIA algorithm.}
\usage{
ilincsSPIA(de = NULL, all = NULL, organism = "hsa",
  data.dir = NULL, pathids = NULL, nB = 2000, plots = FALSE,
  verbose = TRUE, beta = NULL, combine = "fisher")
}
\arguments{
\item{de}{A numeric vector of all differentially expressed (de) genes log2 fold change values named by gene IDs.}

\item{all}{This is a character vector of all used gene IDs in the platform.}

\item{organism}{A character string. One of c("hsa", "mmu", "rno") for Human, Mouse or Rat pathways respectively.}

\item{data.dir}{Path to previously downloaded RData format of different organisms from KEGG database.}

\item{pathids}{.}

\item{nB}{.}

\item{plots}{.}

\item{verbose}{.}

\item{beta}{.}

\item{combine}{.}

\item{authors}{M. Fazel-Najafabadi}
}
\description{
This function allows you to run SPIA method in a very faster way.
}
\examples{
res <- ilincsSPIA(....)
}
\keyword{...}
