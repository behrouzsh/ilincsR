% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{preprocessing}
\alias{preprocessing}
\title{A function for ...}
\usage{
preprocessing(genelist, experiment, path_to_write, up = 4000,
  down = 1000, window.size = 50, org = "Hs", debuging = FALSE)
}
\arguments{
\item{genelist}{A character string containing a list of Entrez gene IDs separated by commas.}

\item{experiment}{The name of selected experiment from portal as a single character string.}

\item{path_to_write}{This parameter specifies the path which user wants to save the results in.}

\item{up, down}{Specifies the maximum number of up and down regulated genes.}

\item{window.size}{...}

\item{test}{This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.}
}
\description{
This function allows you to ...
}
\examples{
res <- batchTregGRS(....)
}
\keyword{...}
