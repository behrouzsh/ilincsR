
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

get_nbr <- function(db=NULL,exp=NULL,geneid=NULL, debuging=FALSE){
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, port = 4040,host="10.165.4.231", password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')}
	sql <- "select AT.Name from ArrayType AT, Array A,Measurement M, MeasurementList ML,Experiment E where E.Measurement_List_ID=ML.ID and ML.Measurement_ID=M.ID and M.Name= A.Name and A.Array_Type_ID=AT.ID and E.Name='";
	Fsql <- paste(sql,exp,"'",sep="")
	rs <- DBI::dbSendQuery(mycon, Fsql)
	AT <- DBI::fetch(rs, n = -1)
	nbr <- NULL	
	for(g in geneid){
                sql <- "select distinct gene2 from GeneDistance.GenePairs where gene1='"
                Fsql<-paste(sql,g,"'",sep="")
		rs <- DBI::dbSendQuery(mycon, Fsql)
		gene2 <- DBI::fetch(rs, n = -1)
		if(nrow(gene2)==0) next
		for(g2 in gene2[,1]){
			sql <- "select distinct GeneID,Symbol,description from GeneDB.GeneInfo where GeneID ='"
			Fsql <- paste(sql,g2,"'",sep="")
			rs <- DBI::dbSendQuery(mycon, Fsql)
			gsd <- DBI::fetch(rs, n = -1)
			if(nrow(gsd)==0) next
			sql <- "select distinct Accession from GeneDB.Genes where Gene='"
			Fsql<-paste(sql,gsd[1,1],"'",sep="")
			rs <- DBI::dbSendQuery(mycon, Fsql)
			acc <- DBI::fetch(rs, n = -1)
			if(nrow(acc)==0) next
			for(a in acc[,1]){
				sql <- "select Gene.ID,Gene.Name from Gene,Feature,Reporter,ArrayType where Gene.Name='"
				sql1 <- "' and Reporter.Gene_List_ID=Gene.ID and Feature.Reporter_ID=Reporter.ID and Feature.Array_Type_ID=ArrayType.ID and ArrayType.Name='"
				Fsql <- paste(sql,a,sql1,AT[1,1],"'",sep="")
				rs <- DBI::dbSendQuery(mycon, Fsql)
				gidname <- DBI::fetch(rs, n = -1)
				if(nrow(gidname)==0) next
				nbr <- c(nbr,g2)
			}
		}

	}
	DBI::dbDisconnect(mycon)
	return(nbr)
}
