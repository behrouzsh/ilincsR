
#' A function to create a modified gct file
#'
#' This function allows you to ...
#' @param esetName Original data in the form of eset.
#' @param path_to_write This parameter specifies the path to save the "gct" format of the downloaded data. 
#'	The default is set to a temp folder: "/mnt/raid/tmp/"
#' @param author M. Fazel-Najafabadi
#' @keywords Download.
#' @export 
#' @examples
#' ## not run
#' mygct <- saveCgMatrixFromEset("some_eset")
#' ## end not run

saveCgMatrixFromEset <- function(esetName,path_to_write="/mnt/raid/tmp/"){
  # path_to_write="/mnt/raid/tmp/"
  # esetName="EDS-1014_GPL10999syn2347004_Wed_Oct__5_15_42_39_2016"
  # esetName="filteredeset_Mon_Aug_7_13_36_48_2017_2313901"
  # esetName="filteredeset_Thu_Aug_10_15_06_05_2017_7218260"
#   library(ArrayTools)
  esetFileName <- paste(path_to_write,esetName,".RData",sep="")
  assign('eset', get(load(esetFileName)))
  eset <- get('eset')
  ex <- Biobase::exprs(eset)
  # pd <-pData(eset)
  # fd <-fData(eset)
  matrixName <- paste(esetName,".txt",sep="")
  matrixPath <- paste(path_to_write,matrixName,sep="")
  # write.matrix(ex, file = matrixPath, sep = "\t")
  # class(colnames(ex))
  # header<-c("rowame",colnames(ex))
  write.table(ex, file = matrixPath, sep = "\t",col.names = NA)
  print(paste(matrixName,sep=""))
}
