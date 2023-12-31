% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigEnrichment.R
\name{sigEnrichment}
\alias{sigEnrichment}
\title{A function for ...}
\usage{
sigEnrichment(sigFile, libName,
  userUploadedProfilePath = "/mnt/raid/tmp/", logPcut = 10,
  sigDB = "ilincs_sigs", two.sided = FALSE, org = "Hs",
  debuging = FALSE)
}
\arguments{
\item{sigFile}{A signature file from ilincs.}

\item{libName}{The library user wants to find enrichments from (currently LIB_5 for drugCategories and LIB_6 for genesCategories).}

\item{userUploadedProfilePath}{Path to the signature file.}

\item{sigDB}{The database for signature libraries ("sigs or new").}

\item{two.sided}{How to calculate p-values for Zscores.}

\item{authors}{M. Fazel-Najafabadi}
}
\description{
This function allows you to ...
}
\examples{
res <- sigEnrichment(sigFile="sig_Fri_Sep_28_16_38_27_2018_2976102.xls", libName="LIB_6", debuging=T)

}
\keyword{mehdi_wtcor}
