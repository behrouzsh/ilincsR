
#' A function for ...
#'
#' This function allows you to ...
#' @param genelist A character string containing a list of Entrez gene IDs separated by commas.
#' @param experiment The name of selected experiment from portal as a single character string.
#' @param path_to_write This parameter specifies the path which user wants to save the results in.
#' @param up,down Specifies the maximum number of up and down regulated genes.
#' @param window.size ...
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

preprocessing <- function(genelist, experiment, path_to_write, up=4000, down=1000, window.size=50, org="Hs", debuging=FALSE)
{
	print(paste("Preprocess called with following parameters at",Sys.time()))
	print(paste("genelist",genelist))
	print(paste("experiment",experiment))
	print(paste("path_to_write",path_to_write))
	print(paste("up",up))
	print(paste("down",down))
	print(paste("window.size",window.size))

##	library(RMySQL)
	requireNamespace("Biobase", quietly = TRUE)
	## First get homologenes
	genes <- find_genes_hgenes_in_platform(exp=experiment, glist=genelist, org=org, debuging=debuging)

	if(is.null(genes) | genes == "")
	{
		return(list(c("No genes found","NA","NA"),NULL))
	}
	#print(genes)
        
	## Get the database name for given experiment first
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
        sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",experiment,"'",sep="")
        rs <- DBI::dbSendQuery(mycon, sql)
        tmp <- DBI::fetch(rs,n=-1)
        DBI::dbDisconnect(mycon)

        db <- tmp[1,"Platform"]
	dataFormat <- tmp[1,"DataFormat"]
        if(is.null(db))
        {
		return(list(c("Experiment not found. Please check the name and try again.","NA","NA"),NULL))
        }
	
	if(!is.null(genes))
	{
		if(dataFormat == "MaxD")
		{
			esetName <- write_eset(genes,experiment,path_to_write, debuging)
		} else {
			esetName <- create_gff_eset(glist=genes,experiment,up,down,window.size,path_to_write, debuging)
		}
	}
	if(esetName == "Experiment not found. Please check the name and try again.")
	{
		return(list(c("Experiment not found. Please check the name and try again.","NA","NA"),NULL))

	}

	if(esetName == "No genes found")
	{
		return(list(c("No genes found","NA","NA"),NULL))

	}
	val <- c(esetName,db,genes)
	#print(val)
	## Get the property values for given experiment
        propVals <- NULL
	load(paste(path_to_write, esetName, ".RData", sep=""))
	eset <- get(esetName)
	for(p in 1:dim(Biobase::pData(eset))[2])
	{
		 propVals <- rbind(propVals,c(colnames(Biobase::pData(eset))[p],paste(unique(Biobase::pData(eset)[,p]),collapse=",",sep="") ))
	}
	colnames(propVals) <- c("Property","Unique values")
	returnList <- list(val,propVals)
	return(returnList)
	
}
