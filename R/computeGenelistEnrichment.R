
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

computeGenelistEnrichment <- function(queryGenelist,libName, ncors=10, debuging=FALSE, sigDB="ilincs_sigs", chost="gimm12.ketl.uc.edu", dhost="gimm2.ketl.uc.edu", dport=3306, org="Hs")
{
	homoloArr <- getHomologousGenes(queryGenelist, org=org)
        homoloArr <- homoloArr[which(!is.na(homoloArr[,3])),]
	if(length(homoloArr) == 2) { # only 1 row
		homoloArr <- t(as.matrix(homoloArr))
	}
	if(dim(homoloArr)[1] == 0){
		return("No homolo genes")
	}
	servSet <- getServerSettings(debuging=debuging, sigDB=sigDB, chost=chost, dhost=dhost, dport=dport)
        mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=servSet$sigDB, host=servSet$host, port = servSet$port, password="public")
        sql <- paste("select pvaluesTablePath,probsTablePath from signatureLibraries where libraryID = '",libName,"'",sep="")
        rs <- DBI::dbSendQuery(mycon, sql)
        dt <- DBI::fetch(rs, n = -1)
        libraryPath <- dt[1,"pvaluesTablePath"]
	libraryScore <- dt[1,"probsTablePath"]
        DBI::dbDisconnect(mycon)
	# load Rdata object for the respective data
	if(libName=="LIB_3") { 
	    refProfiles <- get(load(libraryScore))
	    refgenes <- rownames(refProfiles)
	} else if(libName != "BING") { 
	    refProfiles <- get(load(libraryPath))
	    refgenes <- rownames(refProfiles)
	}
	if(libName=="LIB_3"){ 
	    refData <- -log10(1 - refProfiles)
	} else if(libName != "BING") {
	    refData <- -log10(refProfiles)
	}
# # 	    tt.result <- as.data.frame(data.table::rbindlist(tmp))
	if(libName == "BING") {
	    chopsP <- paste0("cpP_", LETTERS[1:10])
	    tmp <- parallel::mclapply(1:length(chopsP), function(j) {
		refData <- get(chopsP[j], envir=.GlobalEnv)
		tt.result <- apply(refData, 2, function(x) myAllez(x, refgenes, homoloArr[,3]))
		tt.result <- t(tt.result)
		tt.result <- data.frame(dataset=colnames(refData), tt.result, stringsAsFactors=FALSE)
		tt.result
	    }, mc.cores=ncors)
	    tt.result <- as.data.frame(data.table::rbindlist(tmp))
	} else {
	tt.result <- apply(refData,2,function(x) myAllez(x,refgenes,homoloArr[,3]))
	tt.result <- t(tt.result)
	tt.result <- data.frame(dataset=colnames(refData),tt.result, stringsAsFactors=FALSE)
	#return significant profiles and P values
	}
	
	sigResult <- tt.result[ which((tt.result[,3]) > 2.33 ), ,drop=FALSE]
	if(dim(sigResult)[1] == 0){
		sigResult <- tt.result[1:10,]
	}
	sigResult <- sigResult[order(sigResult[,3],decreasing=T),] # if (dim(sigResult)[1] > 0) 
	return(sigResult)
}
