% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadFullSignature.R
\name{loadFullSignature}
\alias{loadFullSignature}
\title{A function to retrieve a previously created full signature}
\usage{
loadFullSignature(sessionID, path_to_write = "/mnt/raid/tmp/")
}
\arguments{
\item{sessionID}{the initial session ID.}

\item{path_to_write}{This optional parameter specifies the path to load the signature from. It should be the same as the path which signature is saved in.}

\item{Value}{The output of this function is in the form of a data.frame. First column is geneIDs, second column is geneNames, 
third column is coefficients and the last column is their probability values for each gene from statistical analysis (ttest).}

\item{authors}{M. Fazel-Najafabadi}
}
\description{
This function allows you to load/retrieve a full signature which is created in the initial step of analysis.
}
\examples{
## not run
fullSig <- loadFullSignature("Fri_Jun_9_15_21_03_2017_4119651")
## end not run
}
\keyword{SignatureTable.}
