
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

random_genes_eset <- function(genes,exp,db,up=4000,down=1000,window.size=100,path_to_write, debuging=FALSE) {
##        library(RMySQL)
	randomGeneEset <- "No genes found"
        geneid <- unlist(strsplit(genes,","))
	noOfGenes <- length(geneid)
        genes <- paste("'",paste(unique(geneid),collapse="','"),"'",sep="")
        sql <- paste("select PG.Gene from ProbeGene PG, Platforms P where P.Database = '",db,"' and P.ID = PG.Platform_ID and PG.Gene not in (",genes,") AND PG.Gene NOT LIKE '%/%' AND PG.Gene IS NOT NULL AND PG.Gene!=''",sep="")
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
        rs <- DBI::dbSendQuery(mycon, sql)
        genesdf <- DBI::fetch(rs, n = -1)
        DBI::dbDisconnect(mycon)
        
        genesdf <- genesdf[!is.na(as.integer(genesdf$Gene)),]
        if(length(genesdf) > 1) {
		randomGenes <- unique(sample(genesdf, noOfGenes, replace=ifelse(length(genesdf) >= noOfGenes, FALSE, TRUE))) ## Mehdi: unique
		randomGenes <- paste(randomGenes,collapse=",",sep="")
		randomGenes <- gsub(",,",",",randomGenes)
		## Get the dataFormat for given experiment first
	# 	servSet <- getServerSettings(test) ## Mehdi: we can remove lines and line 43, 45 & 46
	# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
	# #        if (!test) {
	# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
	# #	else {
	# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
	#         sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
	#         rs <- DBI::dbSendQuery(mycon, sql)
	#         tmp <- DBI::fetch(rs,n=-1)
	#         DBI::dbDisconnect(mycon)
	#         dataFormat<-tmp[1,"DataFormat"]
		dataFormat <- expInfo(exp, debuging)$DataFormat
		if(dataFormat == "MaxD") {
			randomGeneEset <- write_eset(randomGenes,exp,path_to_write, debuging=debuging)
		} else if(dataFormat == "Gff") {
			randomGeneEset <- create_gff_eset(glist=randomGenes,exp,up,down,window.size,path_to_write, debuging=debuging)	
		}
	}
        return(randomGeneEset)
}
