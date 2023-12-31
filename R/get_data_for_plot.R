
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

get_data_for_plot <- function(db,exp,glist, debuging=FALSE)
{
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='kausti', dbname="QueryTables", port = 4040,host="10.165.4.231", password='shinde')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='kausti', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='shinde')}
        sql <- paste("select Platform from ExperimentMetadata where Experiment='",exp , "'",sep="")
        rs <- DBI::dbSendQuery(mycon, sql)
        ar <- DBI::fetch(rs,n=-1)
        arraytp <- ar[1,1]

        sql <- paste("select ProbeGene.Probe as Probe, ProbeGene.Gene as Gene from ProbeGene ,Platforms where ProbeGene.Gene in (",glist,") and Platforms.Platform = '",arraytp,"' and ProbeGene.Platform_ID = Platforms.ID and Platforms.Database ='",db,"'",sep="")
        rs <- dbSendQuery(mycon, sql)
        probegene <- DBI::fetch(rs,n=-1)
        pg <- NULL

        for(i in 1:dim(probegene)[1])
        {
                if(is.null(pg)){
                                pg =  paste(probegene[i,1],":",probegene[i,2],sep="")
                 } else {
                                pg =  paste(pg,",",probegene[i,1],":",probegene[i,2],sep="")
                 }
        }
        DBI::dbDisconnect(mycon)
        return(pg)
}
