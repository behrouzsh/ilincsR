
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

get_gene_probe <- function(db,probegene, debuging=FALSE){
        pg <- unlist(strsplit(probegene,","))
        pg_dataframe <- NULL
        for(p in 1:length(pg)){
                p_g <- unlist(strsplit(pg[p],":"))
                pg_dataframe <- rbind(pg_dataframe,c(p_g[1],p_g[2]))
        }
        f <- get_feature_dataframe(db,pg_dataframe[,1], debuging)
        pg_dataframe <- cbind(pg_dataframe,f[,1])
        colnames(pg_dataframe) <- c("probe","gene","probe_name")
        genes <- unique(pg_dataframe[,2])
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDB", port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", ,host="gimm2.ketl.uc.edu", password='public')}
        gInfo <- NULL
	s <- NULL
	for(a in genes)
	{
		if(is.null(s))
		{
			s <- paste("'",a,"'",sep="")
		} else {
			s <- paste(s,",'",a,"'",sep="")
		}

	}
                sql <- paste("select distinct GeneID,Symbol,description from GeneInfo where GeneID in (",s,")",sep="")
                rs <- DBI::dbSendQuery(mycon, sql)
                geneinfo <- DBI::fetch(rs, n = -1)
        colnames(geneinfo) <- c("gene","symbol","description")

        gp<-merge(pg_dataframe,geneinfo,by.x="gene",by.y="gene")
	DBI::dbDisconnect(mycon)

	return(gp)
}
