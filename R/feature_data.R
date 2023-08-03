
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
#' @keywords test
#' @export 
#' @examples
#' res <- batchTregGRS(....)

feature_data <- function(db,exp,feature,meas,meas_limits, debuging=FALSE){
	## Check if DB was reconstructed
	reconstructed <- FALSE
	dataTable <- NULL
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')
	sql <- paste("select * from ReconfiguredPlatforms where Platform = '",db,"'",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	tmp <- DBI::fetch(rs,n=-1)
	if(nrow(tmp) == 0) reconstructed = FALSE else reconstructed <-  TRUE
	DBI::dbDisconnect(mycon)

    servSet <- getServerSettings(debuging)
    mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db , port = servSet$port, host = servSet$host, password = "public")
	
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')}
	if(reconstructed)
	{
		## Get data table
		sql <- paste("select DataTable from ExperimentDataTable where ExperimentName = '",exp,"'",sep="")
		rs <- DBI::dbSendQuery(mycon, sql)
		dataTable <- DBI::fetch(rs,n=-1)[1,1]
	}
	
        if(!is.null(meas_limits)){
                mnM<-meas_limits[1]
                mxM<-meas_limits[2]
        }
        count<-1
        if(!is.null(meas)){
                ms<-paste(meas,collapse=",")
        }
        s<-paste("'",paste(feature[,2],collapse="','"),"'",sep="")
        if(!is.null(meas)){
		if(!reconstructed) {
                Fsql1<-paste("select D.Value,D.Measurement_ID,F.Name from DoubleProperty D,Feature F where D.Feature_ID in (",s,") and D.Measurement_ID in (",ms,") and D.Feature_ID = F.ID order by F.Name,D.Measurement_ID",sep="")
		} else {
                Fsql1<-paste("select D.Value,D.Measurement_ID,F.Name from ",dataTable," D,Feature F where D.Feature_ID in (",s,") and D.Measurement_ID in (",ms,") and D.Feature_ID = F.ID order by F.Name,D.Measurement_ID",sep="")
		}
        }
        if(!is.null(meas_limits)){
		if(!reconstructed) {
                Fsql1<-paste("select D.Value,D.Measurement_ID,F.Name from DoubleProperty D,Feature F where D.Feature_ID in (",s,") and D.Measurement_ID between ",mnM," and ",mxM," and D.Feature_ID = F.ID order by F.Name,D.Measurement_ID",sep="")
		} else {
                Fsql1<-paste("select D.Value,D.Measurement_ID,F.Name from ",dataTable," D,Feature F where D.Feature_ID in (",s,") and D.Measurement_ID between ",mnM," and ",mxM," and D.Feature_ID = F.ID order by F.Name,D.Measurement_ID",sep="")
                }
        }
        rsFeat <- DBI::dbSendQuery(mycon, Fsql1)
        Feat60 <- DBI::fetch(rsFeat, n = -1)
        colnames(Feat60) <- c("Value","Measurement_ID","Feature_Name")
        unqFeats <- unique(Feat60[,"Feature_Name"])
        nfeature_data <- matrix(Feat60[,"Value"],ncol=length(unqFeats))
        colnames(nfeature_data) <- unqFeats

        DBI::dbDisconnect(mycon)
        return(nfeature_data)
}
