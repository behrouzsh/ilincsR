
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

find_genes_hgenes_in_platform <- function(exp, glist, org="Hs", debuging=FALSE)
{
##	library(RMySQL)
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# 
# #	if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
# 	sql <- paste("select Organism,Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
# 	rs <- DBI::dbSendQuery(mycon, sql)
# 	tmp <- DBI::fetch(rs,n=-1)
# 	DBI::dbDisconnect(mycon)
# 	print(exp)
# 	print(tmp)
# 	dataFormat<-tmp[1,"DataFormat"]
# 	org<-tmp[1,"Organism"]
# 	platform<-tmp[1,"Platform"]
	
	tmp <- expInfo(exp, debuging)
	dataFormat <- tmp$DataFormat
# 	org <- tmp$Organism
	platform <- tmp$Platform
	
	tax <- org2tax(org)
# 	if(org == "human") tax <- "9606"
# 	if(org == "mouse") tax <- "10090"
# 	if(org == "rat") tax <- "10116"
		
	glist <- unlist(strsplit(glist,","))
	glist <- sub(" ","",glist)
	glist <- gsub(",,",",",glist)
        glist <- gsub("'NA',","",glist)
        glist <- gsub(",'NA'","",glist)

# 	load(url("http://www.eh3.uc.edu/geneinfo.RData"))
	data(geneinfo)
	genedf <- geneinfo[match(glist,geneinfo[,1]),]
	tgenedf <- genedf[genedf["TaxID"]==tax,]

	## homologenes
	hml <- genedf[genedf["TaxID"]!=tax,"HomologeneID"]
	x <- geneinfo[geneinfo["TaxID"]==tax,]
	ind <- match(hml,x[,"HomologeneID"])
	hgenedf <- x[ind,]

	## genes for experiment's organism
	finalgenes <- c(hgenedf[,"GeneID"],tgenedf[,"GeneID"])
	finalgenes <- paste("'",paste(finalgenes,collapse="','"),"'",sep="")
	finalgenes <- gsub(",,",",",finalgenes)

        finalgenes <- gsub("'NA',","",finalgenes)
        finalgenes <- gsub(",'NA'","",finalgenes)

	## genes in the platform
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")

#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
	sql <- paste("SELECT distinct Gene FROM ProbeGene PG,Platforms P where P.Platform ='",platform,"' AND PG.Platform_ID = P.ID AND Gene in (",finalgenes,") AND Gene NOT LIKE \"%/%\"",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	tmp <- DBI::fetch(rs,n=-1)
	platformgenes <- paste(tmp[!is.na(as.integer(tmp[,1])), 1],collapse=",",sep="")
	platformgenes <- gsub(",,",",",platformgenes)
	DBI::dbDisconnect(mycon)
	return(platformgenes)	
}
