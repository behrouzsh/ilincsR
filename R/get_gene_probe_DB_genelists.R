
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

get_gene_probe_DB_genelists <- function(db,genelistid,platform, debuging=FALSE){
##	library(RMySQL)
	servSet <- getServerSettings(debuging)
	GDBcon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDB", port = servSet$port, host = servSet$host, password = "public")
	GDistcon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDistance", port = servSet$port, host = servSet$host, password = "public")
	dbcon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  GDBcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", port = 4040,host="10.165.4.231", password='public')
#	  GDistcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDistance", port = 4040,host="10.165.4.231", password='public')
#	  dbcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  GDBcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", host="gimm2.ketl.uc.edu", password='public')
#	  GDistcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDistance", host="gimm2.ketl.uc.edu", password='public')
#	  dbcon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')}
	## get platform ID
	sql <- paste("select ID from ArrayType where Name='",platform,"'",sep="")
	rs <- DBI::dbSendQuery(dbcon, sql)
	platformid <- DBI::fetch(rs, n = -1)

	sql <- paste("select distinct ProbeName,ProbeID,GeneID from GeneListsProbes where PlatformID =",platformid[1,1]," and GeneListID=",genelistid) 
	rs <- DBI::dbSendQuery(GDistcon, sql)
        geneprobe <- DBI::fetch(rs, n = -1)
	colnames(geneprobe) <- c("probe_name","probe","gene")
	gInfo <- NULL
        for(g in 1:dim(geneprobe)[1]){
                sql <- paste("select distinct GeneID,Symbol,description from GeneInfo where GeneID = '",geneprobe[g,3],"'",sep="")
                rs <- DBI::dbSendQuery(GDBcon, sql)
                geneinfo <- DBI::fetch(rs, n = -1)
                if(nrow(geneinfo)==0) next
                gInfo <- rbind(gInfo,c(geneinfo[1,1],geneinfo[1,2],geneinfo[1,3]))
        }
        colnames(gInfo) <- c("gene","symbol","description")

        gp <- merge(geneprobe,gInfo,by.x="gene",by.y="gene")
	
	DBI::dbDisconnect(dbcon)
	DBI::dbDisconnect(GDBcon)
	DBI::dbDisconnect(GDistcon)

        return(gp)
}
