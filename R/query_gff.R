
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

query_gff <- function(experiment, genelist, up, down, debuging=FALSE, smpls="all")
{
##	library(RMySQL)

	## IMP: We are deliberately connecting to a specific database instead of load balancer here because we are creating views and subsequent queries need to use views on the same server
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# #@!	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname='QueryTables', host='gimm2.ketl.uc.edu', password='public')
# #@	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname='QueryTables', host='db2.ketl.uc.edu', password='public')
# 	sql <- paste("select Platform from ExperimentMetadata where Experiment='",experiment,"'",sep="")
# # 	rs <- DBI::dbSendQuery(mycon, sql)
# # 	db <- DBI::fetch(rs,n=-1)[1,1]
# 	db <- DBI::dbGetQuery(mycon, sql)[1,1]
# 	DBI::dbDisconnect(mycon)
	db <- expInfo(experiment, debuging)$Platform

	# Code changed by vkj on : Fri Oct 21 16:21:32 EDT 2011
	# mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host='db2.ketl.uc.edu', password='public')
##	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='ilincs_tmp', dbname=db, host='db5.ketl.uc.edu', password='shinde')
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='ilincs_tmp', dbname=db, host='gimm2.ketl.uc.edu', password='shinde')
##	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='kausti', dbname=db, host='gimm2.ketl.uc.edu', password='shinde')
#@@	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='kausti', dbname=db, host='db3.ketl.uc.edu', password='shinde')
#@	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public2', dbname=db, host='db2.ketl.uc.edu', password='public2')

	sessionID <- paste(c(strsplit(date(),split=" ")[[1]],as.integer(runif(1)*10e6)),sep="",collapse="_")
	sessionID <- gsub(":","_",sessionID)

	##Create view given gene list
	sql <- paste("create or replace view tmp.",sessionID,
		    " as SELECT RefID,GeneID,Chromosome,Strand,txStart,txEnd,txStart as TSS FROM refGene WHERE  Strand = '+' AND GeneID in   (", 
		    genelist,") Union SELECT RefID,GeneID,Chromosome,Strand,txStart,txEnd,txEnd as TSS FROM refGene WHERE  Strand = '-' AND GeneID in (", 
		    genelist,") Group by Chromosome order by TSS", sep="")

	rs <- DBI::dbSendQuery(mycon, sql)

	##Get data for TSS +- range bps

	## Get samples for the experiment
	sql <- paste("select S.Sample_Name,DS.Value_Column from Dataset D, Sample S,Dataset_Sample DS where D.Dataset_Name = '",experiment,
			"' and D.Dataset_ID = DS.Dataset_ID and DS.Sample_ID = S.Sample_ID",sep="")
# 	rs <- DBI::dbSendQuery(mycon, sql)
# 	samples <- DBI::fetch(rs,n=-1)
	samples <- DBI::dbGetQuery(mycon, sql)
#	samples <- samples
	if (smpls != "all") {
# 	  smpls <- unlist(strsplit(smpls, ","))
	  samples <- samples[(samples[, "Sample_Name"] %in% unlist(strsplit(smpls, ","))), ,drop=FALSE]
	}
	
	resultdf <- NULL
	for(sample in 1:length(samples[,"Sample_Name"]))
	{
		print(paste("sample",samples[sample, 2]))
		sql <- paste("SELECT DISTINCT t.*,hg.chrom,hg.chromStart,hg.chromEnd,hg.",samples[sample,2], 
			    ",abs(t.TSS-hg.chromStart) as startDiff,abs(t.TSS-hg.chromEnd) as endDiff  FROM tmp.",sessionID, 
			    " t,",samples[sample,1],
			    " hg WHERE hg.chrom = t.Chromosome AND t.Strand = \"+\" AND ( hg.chromStart BETWEEN t.TSS - ",up, 
			    " AND t.TSS + ", down, " OR hg.chromEnd BETWEEN t.TSS - ", up, " AND t.TSS + " , down, 
			    " OR hg.chromStart <= t.TSS AND hg.chromEnd >= t.TSS ) Union SELECT DISTINCT t.*,hg.chrom,hg.chromStart,hg.chromEnd,hg.", 
			    samples[sample,2],",abs(t.TSS-hg.chromStart) as startDiff,abs(t.TSS-hg.chromEnd) as endDiff  FROM tmp.",
			    sessionID," t, ",samples[sample,1] , 
			    " hg WHERE hg.chrom = t.Chromosome AND t.Strand = \"-\" AND ( hg.chromStart BETWEEN t.TSS - ",down, 
			    " AND t.TSS + ", up, " OR hg.chromEnd BETWEEN t.TSS - ", down, " AND t.TSS + " , up, 
			    " OR hg.chromStart <= t.TSS AND hg.chromEnd >= t.TSS )",sep="")
#		print(sql)
# 		rs <- DBI::dbSendQuery(mycon, sql)
# 		returndf <- DBI::fetch(rs,n=-1)
		returndf <- DBI::dbGetQuery(mycon, sql)
		if(nrow(returndf)!=0){
			returndf <- cbind(returndf,rep(samples[sample,1],dim(returndf)[1]))
			#colnames(returndf) <- c(colnames(returndf),samples[sample,"Sample_Name"])
			resultdf <- rbind(resultdf,returndf)
		}
	}
	sql <- paste("drop view tmp.",sessionID,sep="")
	rs <- DBI::dbSendQuery(mycon, sql)

	DBI::dbDisconnect(mycon)
	if(!is.null(resultdf)){
		colnames(resultdf)[dim(resultdf)[2]] <- "Sample"
	}
	return(resultdf)

}
