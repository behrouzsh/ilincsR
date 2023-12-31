% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigPlot.R
\name{sigPlot}
\alias{sigPlot}
\title{A function for ...}
\usage{
sigPlot(sigFile = NULL, sigID = NULL,
  path_to_write = "/mnt/raid/tmp/", output = "json", dotSize = 3,
  fontSize = 15, debuging = FALSE)
}
\arguments{
\item{sigFile}{This should be a signature file created or uploaded by user.}

\item{sigID}{This is a precomputed signature from any sig-library.}

\item{path_to_write}{This parameter specifies the path which user wants to save the results.}

\item{debuging}{This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
alternative database server "gimm2" if debuging=TRUE dev servers will be used.}

\item{display}{Display plot or just get the object.}

\item{authors}{M. Fazel-Najafabadi}
}
\description{
This function allows you to ...
}
\examples{
## Do not run
sigFile <- "sig_Thu_Apr_12_12_45_57_2018_298226.xls"
sigID <- "LINCSCP_109"
path_to_write <- "/opt/raid10/genomics/mehdi/ilincs/output/"
plt <- sigPlot(sigFile, sigID, path_to_write, debuging=TRUE, output=TRUE)
## End do not run
}
\keyword{filter}
