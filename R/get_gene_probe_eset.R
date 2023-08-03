
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

get_gene_probe_eset <- function(db, probes, debuging=FALSE)
{
	probeNames <- NULL
	probeNames <- paste("'", paste(probes, collapse="','", sep=""), "'", sep="")

	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
	sql <- paste("select PG.Probe,PG.Gene from ProbeGene PG,Platforms PL where PL.Platform ='",db,"' and PG.Platform_ID = PL.ID and PG.Probe in (",probeNames,")",sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	pg_dataframe <- DBI::fetch(rs,n=-1)
	colnames(pg_dataframe) <- c("probe_name","gene")
	pg_dataframe <- pg_dataframe[!is.na(as.integer(pg_dataframe$gene)), ] ##Dec,2017
	organism <- DBI::dbGetQuery(mycon, paste0("select Organism from ExperimentMetadata where Platform='", db, "'"))[1,1]
	tax <- switch(organism, human = 9606, mouse = 10090, rat = 10116)
	DBI::dbDisconnect(mycon)

        genes <- unique(pg_dataframe[,"gene"])
	s <- NULL
	s <- paste("'",paste(genes,collapse="','",sep=""),"'",sep="")

	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDB", port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", host="gimm2.ketl.uc.edu", password='public')}
	sql <- paste("select distinct GeneID, Symbol, description from GeneInfo where GeneID in (",s,") and TaxID='", tax, "'", sep="")
	rs <- DBI::dbSendQuery(mycon, sql)
	geneinfo <- DBI::fetch(rs, n = -1)
        colnames(geneinfo) <- c("gene","symbol","description")
# 	geneinfo <- geneinfo[!is.na(geneinfo$gene), ]
	DBI::dbDisconnect(mycon)

# 	geneinfo <- geneid2symbol(paste0(genes, collapse=","))
# 	colnames(geneinfo) <- c("gene","symbol","description")
	
	gp <- merge(pg_dataframe, geneinfo, by.x="gene", by.y="gene")
	rownames(gp) <- gp$probe_name

        return(gp)
}
