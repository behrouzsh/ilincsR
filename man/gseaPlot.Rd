% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gseaPlot.R
\name{gseaPlot}
\alias{gseaPlot}
\title{A function for ...}
\usage{
gseaPlot(sigFile = NULL, sigID = NULL, glist,
  path_to_write = "/mnt/raid/tmp/", output = c("png", "pdf"),
  fontSize = 15, debuging = FALSE)
}
\arguments{
\item{sigID}{This is a precomputed signature from any sig-library.}

\item{glist}{A list of Entrez gene IDs to enrich.}

\item{path_to_write}{This parameter specifies the path which user wants to save the results.}

\item{output}{Type of the output plot including "png", "pdf", "ggplot", "json" and "widget".}

\item{debuging}{This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
alternative database server "gimm2" if debuging=TRUE dev servers will be used.}

\item{sigfile}{This should be a signature file created or uploaded by user.}

\item{authors}{M. Fazel-Najafabadi}
}
\description{
This function allows you to ...
}
\examples{
## Do not run
sigfile <- "sig_Thu_Apr_12_12_45_57_2018_298226.xls"
sigID <- "LINCSCP_109"
path_to_write <- "/opt/raid10/genomics/mehdi/ilincs/output/"
plt <- sigPlot(sigfile, sigID, path_to_write, debuging=TRUE, output=TRUE)
## End do not run
}
\keyword{filter}
