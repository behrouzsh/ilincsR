
#' A function to retrieve an experiment metadata
#'
#' This function allows you to retrieve necessary metadata about an experiment from databases.
#' @param exp The experiment which you want to get metadata for.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param author M. Fazel-Najafabadi
#' @keywords Download.
#' @export 
#' @examples
#' ## not run
#' EDS1014 <- expInfo(exp="EDS-1014")
#' ## end not run

expInfo <- function(exp, debuging=FALSE) {
    servSet <- getServerSettings(debuging)
    mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
    expinfo <- DBI::dbGetQuery(mycon, paste0("select * from ExperimentMetadata where experiment = '", exp, "'"))
    DBI::dbDisconnect(mycon)
    
    Platform <- expinfo$Platform
    dataTable <- expinfo$dataTable
    tblName <- expinfo$tblName
    dtRows <- expinfo$dtRows
    lincs_dsgc <- expinfo$lincs_dsgc
    Organism <- expinfo$Organism
    SampleType <- expinfo$SampleType
    DataType <- expinfo$DataType
    Assay <- expinfo$Assay
    Description <- expinfo$Description
    PubMedDescription <- expinfo$PubMedDescription
    GeoLink <- expinfo$GeoLink
    PubLink <- expinfo$PubLink
    DataFormat <- expinfo$DataFormat
    SourceID <- expinfo$SourceID
    eset_path <- expinfo$eset_path
    Rprogram_path <- expinfo$Rprogram_path
    analyst <- expinfo$analyst
    
    info <- list(Platform=Platform, dataTable=dataTable, tblName=tblName, dtRows=dtRows, lincs_dsgc=lincs_dsgc, Organism=Organism, 
		 SampleType=SampleType, DataType=DataType, Assay=Assay, Description=Description, PubMedDescription=PubMedDescription, 
		 GeoLink=GeoLink, PubLink=PubLink, DataFormat=DataFormat, SourceID=SourceID, eset_path=eset_path, Rprogram_path=Rprogram_path, analyst=analyst
		  )
    return(info)
}

