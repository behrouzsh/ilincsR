% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloadData.R
\name{downloadData}
\alias{downloadData}
\title{A function to download whole or part of a dataset from ilincs.org portal}
\usage{
downloadData(exp, filterchk = NULL, includeORexclude = "1",
  glist = NULL, path_to_write = "/mnt/raid/tmp/", display = TRUE,
  homoloGenes = TRUE, id.na.rm = TRUE, debuging = FALSE)
}
\arguments{
\item{exp}{The experiment which you want to filter and download. It should be present in the portal's database.}

\item{filterchk}{filterchk is in the form of a string. It can have multiple property of the pData to filter
based on. Each pair of property:value should be separated by a comma ",". Pairs are saparated
by colon ":". An example ot filterchk is :
filterchk="property1:value1,property1:value2,property2:value1"}

\item{includeORexclude}{This argument basically is designed to filter the pData based on what 
filterchk is or selected property of eset when is set to "1", or the reverse selection of filterchk 
when set to "2". 
It also can have NULL or "n" values when there is no property selected, in this case 
the function will return the original ExpressionSet.}

\item{path_to_write}{This parameter specifies the path to save the "gct" format of the downloaded data. 
The default is set to a temp folder: "/mnt/raid/tmp/"}

\item{display}{If user wants to get the gctFile as an R object (text) otherwise it will only be saved and file name will be returned.}

\item{test}{This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.}

\item{author}{M. Fazel-Najafabadi}
}
\description{
This function allows you to filter and download a dataset as "gct" format from ilincs portal.
}
\examples{
## not run
mygct <- downloaData(exp="EDS-1014", filterchk="subtype:Luminal,subtype:Basal", path_to_write="your/path/to/write/here", display=FALSE)
## end not run
}
\keyword{Download.}
