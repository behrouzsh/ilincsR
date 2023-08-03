
#' A function for ...
#'
#' This function allows you to ...
#' @param query.p ...
#' @param reference.p ...
#' @param query.gL ...
#' @param reference.gL ...
#' @param na.rm ...
#' @param estimateNullDistr ...
#' @param nullDistrQuantiles ...
#' @param nullDistrN ...
#' @param tolerateWarnings ...
#' @param pAdjust.method.query ...
#' @param pAdjust.method.reference ...
#' @param lambda ...
#' @param rescale ...
#' @param plotRescaled ...
#' @param multicoreN ...
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

master_for_servlet <- function(db,exp,prop,filterchk,filters,ifCluster,includeORexclude,esetName, 
				ifLR,path_to_write,up,down,window.size, debuging=FALSE)
{
##	source("http://eh3.uc.edu/r/gimmHeat.R")
	load(paste(path_to_write,esetName,".RData",sep=""))
	
	## Get dataType of the experiment
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#@!	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", port = 4040,host="10.165.4.231", password="public")}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", ,host="gimm2.ketl.uc.edu", password="public")}
        sql <- paste("select DataType,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#@!        sql <- paste("select DataType,DataFormat from ExperimentMetadata where experimentName = '",exp,"'",sep="")
        rs <- DBI::dbSendQuery(mycon, sql)
        dt <- DBI::fetch(rs, n = -1)
        dataType <- dt[1,"DataType"]
	dataFormat <- dt[1,"DataFormat"]
        DBI::dbDisconnect(mycon)
	esetString <- esetName
	esetName <- get(esetName)

	meas_limits <- initiate(db,exp,web=TRUE,dataFormat, debuging=debuging)

	

	totalSamples <- NULL
	filteredSamples <- NULL
		
	## If property is specified
	## Mehdi: one line for Gff clustering
	if (dataFormat == "Gff") ifCluster <- "row"
	totalSamples <- dim(Biobase::exprs(esetName))[2]
	probes <- rownames(Biobase::exprs(esetName))
	geneprobes <- get_gene_probe_eset(db, probes, debuging)
	resultdf <- plot_with_properties_eset(esetName, geneprobes, prop,web=FALSE, dataType, filterchk, filters, includeORexclude,
					    ifCluster, ifLR, exp, db, path_to_write, dataFormat, up, down, window.size)
	
	if(resultdf[1,"Remark"] != "Filterd zero" & resultdf[1,"Remark"] != "all missing")
	{
		samples <- resultdf[1,"FilteredSamples"]
		newEset <- write_fitered_eset(esetString, samples, db, exp, path_to_write)	
		resultdf <- cbind(resultdf, newEset)
	}
	return(resultdf)
}
