
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

nbr_struct <- function(uid,glid,db, debuging=FALSE)
{
##	library(RMySQL)
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
	GDistcon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDistance", port = servSet$port, host = servSet$host, password = "public")
	authcon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "authentication", port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, port = 4040,host="10.165.4.231", password='public')
#	  GDistcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDistance", port = 4040,host="10.165.4.231", password='public')
#	  authcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="authentication", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')
#	  GDistcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDistance", host="gimm2.ketl.uc.edu", password='public')
#	  authcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="authentication", host="gimm2.ketl.uc.edu", password='public')}
	expts <- NULL
	df <- NULL

	sql <- paste("select DP.Dataset from User_Permissions upm,Permissions pm, Datasets_Projects DP where upm.User =",uid, " and ( pm.PermissionID = upm.Permissions or pm.PermissionID = 1) and DP.ID = pm.DatasetID order by DP.ID")
	rs <- DBI::dbSendQuery(authcon, sql)
        datasets <- DBI::fetch(rs, n = -1)
	colnames(datasets) <- "datasets"

	sql <- paste("select List from geneLists where ID=",glid)
	rs <- DBI::dbSendQuery(GDistcon, sql)
	list <- DBI::fetch(rs, n = -1)
	glst <- unlist(strsplit(list[1,1],":"))
		
	i <- NULL
	for(i in 1:(length(glst)-1))
	{
		for(j in (i+1):length(glst))
		{
			sql <- paste("select * from GenePairs where gene1 ='",glst[i],"' and gene2 = '",glst[j],"'",sep="") 
			rs <- DBI::dbSendQuery(GDistcon, sql)
			gpairs <- DBI::fetch(rs, n = -1)
			if(nrow(gpairs) == 0) next
			
			expts <- unlist(strsplit(gpairs[1,3],","))
			if(length(expts) == 0) next
			df <- rbind(df,c(gpairs[1,1],gpairs[1,2],length(expts),gpairs[1,3]))
		}
	}
	DBI::dbDisconnect(GDistcon)
	DBI::dbDisconnect(authcon)
	DBI::dbDisconnect(mycon)
	return(df)
}
