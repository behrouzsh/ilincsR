
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

# initiate <- function(db=NULL,exp=NULL,web=TRUE,dataFormat, debuging=FALSE){
# ##        library(DBI)
# ##        library(RMySQL)
# ##        library(ctc,verbose=F)
#         #library(R2HTML,verbose=F)
# ##        library(gplots,verbose=F,warn.conflicts=F)
# ##        library(marray,verbose=F,warn.conflicts=F)
#         pal <- marray::maPalette(low="blue",high="yellow",mid="black")
#         allColors <- c("blue","red","green","yellow","purple","brown","black","blue4","red4","green4","yellow4")
# 
# 	if(dataFormat != "Gff")
# 	{
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
# #        if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=db, host="10.165.4.231", port = 4040, password="public")}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=db, host="gimm2.ketl.uc.edu", password="public")}
# 	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
# 
#         Fsql1 <- paste("select Measurement_List_ID from Experiment where Name='",exp,"'",sep="")
#         rsml <- DBI::dbSendQuery(mycon, Fsql1)
#         mlid <- DBI::fetch(rsml, n = -1)
# 
#         Fsql1 <- paste("select min(Measurement_ID) from MeasurementList where ID=",mlid[1,1])
#         rsmnM <- DBI::dbSendQuery(mycon, Fsql1)
#         mnM <- DBI::fetch(rsmnM, n = -1)
# 
#         Fsql1 <- paste("select max(Measurement_ID) from MeasurementList where ID=",mlid[1,1])
#         rsmxM <- DBI::dbSendQuery(mycon, Fsql1)
#         mxM <- DBI::fetch(rsmxM, n = -1)
# 
#         meas_limits<-c(mnM,mxM)
#         DBI::dbDisconnect(mycon)
#         return(meas_limits)
# 	}
# 	return(NULL)
# }

initiate <- function(db=NULL,exp=NULL,web=TRUE,dataFormat, debuging=FALSE){

	if(dataFormat != "Gff") {
	  servSet <- getServerSettings(debuging)
	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
	  mnM <- DBI::dbGetQuery(mycon, paste0("select min(Measurement_ID) from MeasurementMetadata where Experiment = '", exp, "'"))
	  mxM <- DBI::dbGetQuery(mycon, paste0("select max(Measurement_ID) from MeasurementMetadata where Experiment = '", exp, "'"))
	  meas_limits <- c(mnM,mxM)
	  DBI::dbDisconnect(mycon)
	  return(meas_limits)
	} else {
	  return(NULL)
	}
}
