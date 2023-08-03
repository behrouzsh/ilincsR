
#' A function to download whole or part of a dataset from ilincs.org portal
#'
#' This function allows you to filter and download a dataset as "gct" format from ilincs portal.
#' @param exp The experiment which you want to filter and download. It should be present in the portal's database.
#' @param filterchk filterchk is in the form of a string. It can have multiple property of the pData to filter
#'	based on. Each pair of property:value should be separated by a comma ",". Pairs are saparated
#'	by colon ":". An example ot filterchk is :
#'	filterchk="property1:value1,property1:value2,property2:value1"
#' @param includeORexclude This argument basically is designed to filter the pData based on what 
#'	filterchk is or selected property of eset when is set to "1", or the reverse selection of filterchk 
#'	when set to "2". 
#'	It also can have NULL or "n" values when there is no property selected, in this case 
#'	the function will return the original ExpressionSet.
#' @param path_to_write This parameter specifies the path to save the "gct" format of the downloaded data. 
#'	The default is set to a temp folder: "/mnt/raid/tmp/"
#' @param display If user wants to get the gctFile as an R object (text) otherwise it will only be saved and file name will be returned.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param author M. Fazel-Najafabadi
#' @keywords Download.
#' @export 
#' @examples
#' ## not run
#' mygct <- downloaData(exp="EDS-1014", filterchk="subtype:Luminal,subtype:Basal", path_to_write="your/path/to/write/here", display=FALSE)
#' ## end not run

downloadData <- function(exp, filterchk=NULL, includeORexclude="1", glist=NULL, path_to_write="/mnt/raid/tmp/", display=TRUE, homoloGenes=TRUE, id.na.rm=TRUE, debuging=FALSE) {
#  library(RMySQL)
#  library(Biobase)
#  library(ilincsR)

fileName <- paste(c(exp, strsplit(date(),split=" ")[[1]], as.integer(runif(1)*10e6)),sep="",collapse="_")
# if(is.null(filterchk) & is.null(glist)) fileName <- exp
if(!is.null(filterchk)) if (filterchk=="n") filterchk <- NULL

if(is.null(filterchk) & is.null(glist)) {if(id.na.rm) fileName <- exp else fileName <- paste0(exp, "_with_na")}
    
# fileName <- paste0(gsub(":","_", fileName), ".gct")
fileName <- gsub(":","_", fileName)
res <- list()
res$fileName <- fileName

# servSet <- getServerSettings(test) ## Mehdi: remember to change ir back to test after database sync
# mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# esetFileName <- DBI::dbGetQuery(mycon, paste("select eset_path from ExperimentMetadata where experiment = '", exp, "'", sep=""))
# #DBI::dbClearResult(RMySQL::dbListResults(mycon)[[1]])
# DBI::dbDisconnect(mycon)
#   
# esetFileName <- esetFileName$eset_path[1]
##esetFileName <- expInfo(exp, test)$eset_path
##esetFileName <- gsub("/mnt/", "/opt/", esetFileName)
#   assign(exp, get(load(esetFileName)))
#   eset <- get(exp)
##  if(is.null(filterchk)) eset <- getAnnotatedExp(exp, test) else eset <- get(load(esetFileName))

if(!file.exists(paste0(path_to_write, fileName, ".RData"))) {# paste0(path_to_write, res$fileName)))

	if(is.null(glist)) {
	    eset <- getAnnotatedExp(exp, id.na.rm, debuging)
	    eset <- filter_eset(eset, filterchk=filterchk, includeORexclude=includeORexclude)
	} else {
	    glist <- gsub(",NA","", glist)
	    glist <- gsub("NA,","", glist)
	    eset <- getData(exp, glist=glist, filterchk=filterchk, includeORexclude=includeORexclude, homoloGenes=homoloGenes, debuging=debuging)
	}
  if(dim(Biobase::fData(eset))[2] == 0) {
      fdata <- data.frame(Gene_id=rownames(Biobase::exprs(eset)), row.names=rownames(exprs(eset)))
      Biobase::featureData(eset) <- as(fdata, "AnnotatedDataFrame")
  }
  res$gct <- eset
  
  save(eset, file=paste(path_to_write, ifelse((fileName==exp | fileName==paste0(exp, "_with_na")), "/", "/filteredeset_"), fileName, ".RData", sep=""))
  esetToGct1.3(eset=eset, gctFile=res$fileName, path_to_write=path_to_write)

  if (display) return(res) else return(res$fileName)
} else {
  message("File already exists!")
  return(res$fileName)
}
}



    
    
    
    
    
    
    
    