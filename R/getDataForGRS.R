
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

getDataForGRS <- function(exp,filterchk,filters,includeORexclude,samples,platform, debuging=FALSE)
{
##	library(RMySQL)
	a <- Sys.time()
	## Check if DB was reconstructed
	reconstructed <- FALSE
	dataTable <- NULL
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="10.165.4.231", port = 4040, password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
	sql <- paste("select * from ReconfiguredPlatforms where Platform = '",platform,"'",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	tmp <- DBI::fetch(rs,n=-1)
	if(nrow(tmp) == 0) {
		reconstructed <- FALSE
	}
	else { 
		reconstructed <-  TRUE 
	}
	DBI::dbDisconnect(mycon)

	ms <- paste(paste("'",samples,collapse="',",sep=""),"'",sep="")

	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = platform, port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="10.165.4.231", port = 4040, password="public")}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="gimm2.ketl.uc.edu", password="public")}
	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
	if(reconstructed)
	{
		## Get data table
		sql <- paste("select DataTable from ExperimentDataTable where ExperimentName = '",exp,"'",sep="")
		rs <- DBI::dbSendQuery(mycon, sql)
		dataTable <- DBI::fetch(rs,n=-1)[1,1]
	}
	# Get meas ids
	sql <- paste("Select M.ID,M.Name from Measurement M,MeasurementList ML,Experiment E where E.Name = '",exp,"' AND E.Measurement_List_ID=ML.ID AND ML.Measurement_ID = M.ID AND M.Name in (",ms,")",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	meas <- DBI::fetch(rs,n=-1)
	ms <- paste(meas[,"ID"],collapse=",",sep="")
	if(!reconstructed)
		sql <- paste("select * from DoubleProperty where Measurement_ID in(",ms,")",sep="")
	else
		sql <- paste("select * from ",dataTable," where Measurement_ID in(",ms,")",sep="")	
	c <- Sys.time()
	rs <- DBI::dbSendQuery(mycon, sql)
	d <- Sys.time()
	dp <- DBI::fetch(rs,n=-1)
	e <- Sys.time()
	
	feats <- paste(dp[,"Feature_ID"],collapse=",",sep="")
	sql <- paste("Select * from Feature where ID in (",feats,")",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	ft <- DBI::fetch(rs,n=-1)
	
	o <- match(dp[,"Feature_ID"],ft[,"ID"])
	dp <- cbind(dp[,c("Measurement_ID","Value")],ft[o,"Name"])
	colnames(dp) <- c("Measurement_ID","Value","Probe")
	DBI::dbDisconnect(mycon)
	
	# Get probegene
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="10.165.4.231", port = 4040, password="public")}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="gimm2.ketl.uc.edu", password="public")}
	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="192.168.1.201", port = 3306, password="public")
	sql <- paste("select PG.Probe,PG.Gene from ProbeGene PG, Platforms PT where PT.Platform = '",platform,"' AND PT.ID = PG.Platform_ID",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	pg <- DBI::fetch(rs,n=-1)
	
	DBI::dbDisconnect(mycon)
	
	o <- match(dp[,"Probe"],pg[,"Probe"])
	dp <- cbind(dp,pg[o,"Gene"])
	colnames(dp) <- c("Measurement_ID","Value","Probe","Gene")
	
	o <- match(dp[,"Measurement_ID"],meas[,"ID"])
	dp <- cbind(meas[o,"Name"],dp[,c("Probe","Gene","Value")])
	colnames(dp) <- c("Sample","Probe","Gene","Value")


	r <- dp
	b <- Sys.time()
	print(paste("start time query",a))
	print(paste("end time query",b))
	print(paste("time for query",b-a))

	return(r)
}
