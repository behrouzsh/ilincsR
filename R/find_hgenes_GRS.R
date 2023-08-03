
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

find_hgenes_GRS <- function(exp,org,glist, debuging=FALSE)
{
	if(!is.null(exp))
	{
##		library(RMySQL)
		servSet <- getServerSettings(debuging)
		mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#		if (!test) {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="10.165.4.231", port = 4040, password='public')}
#		else {
#		mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
		#mycon <- dbConnect(MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
		sql <- paste("select Organism,Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
		rs <- DBI::dbSendQuery(mycon, sql)
		tmp <- DBI::fetch(rs,n=-1)
		DBI::dbDisconnect(mycon)
		dataFormat <- tmp[1,"DataFormat"]
		org <- tmp[1,"Organism"]
		platform <- tmp[1,"Platform"]
	}
	tax <- "9606"
	if(org == "human") tax<-"9606"
	if(org == "mouse") tax<-"10090"
	if(org == "rat") tax<-"10116"
		
	glist <- unlist(strsplit(glist,","))
	glist <- sub(" ","",glist)
	glist <- gsub(",,",",",glist)
        glist <- gsub("'NA',","",glist)
        glist <- gsub(",'NA'","",glist)

	load(url("http://www.eh3.uc.edu/geneinfo.RData"))
	genedf <- geneinfo[match(glist,geneinfo[,1]),]
#	tgenedf<-genedf[genedf["TaxID"]==tax,]

	l <- which(genedf[,"TaxID"] == tax)

	## homologenes
	hml <- genedf[genedf["TaxID"]!=tax,"HomologeneID"]
	x <- geneinfo[geneinfo["TaxID"]==tax,]
	ind <- match(hml,x[,"HomologeneID"])
	hgenedf <- x[ind,]
	
	if(length(l) > 0)
	{
		hgenedf[l,] <- genedf[l,]
	}

	retArr <- NULL
	retArr <- cbind(genedf[,"GeneID"], hgenedf[,"GeneID"])
	colnames(retArr) <- c("geneID","homoloGeneID")
	return(retArr)
}
