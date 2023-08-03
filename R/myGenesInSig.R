
#' A function for finding common genes between a list and an experiment.
#'
#' This function allows you to know how many of the genes in your list are or have homologs in a spesified experiment.
#' @param exp ...
#' @param glist ...
#' @param test ...
#' @param authors M. Fazel-Najafabadi
#' @keywords test
#' @export 
#' @examples
#' ## Not run
#' genelist <- "1002,9447,10001,10096,79575"
#' res <- myGenesInExp(exp="EDS-1014", glist=genelist, debuging=FALSE)
#' ## End, not run

myGenesInSig <- function(sigID, glist, org="Hs", sigDB="ilincs_sigs", debuging=FALSE) {

	servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
	glist <- paste(glist, collapse="','")
	
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = servSet$sigDB, port = servSet$port, host = servSet$host, password = "public")
	
	sigType <- strsplit(sigID, split="_")[[1]][1]
	sigTable <- switch(sigType,
#		      LINCSKD = "LincsKdSigData_Nov12",
		      LINCSKD = "LincsKdSigData",
#		      LINCSCP = "LincsCpSigData_Nov12",
		      LINCSCP = "LincsCpSigData",
		      LINCSSH = "LincsShSigData",
		      `EDS-1014` = "LincsRnaSeqSigData",
		      `LDS-1237` = "LincsRnaSeqSigData",
		      `LDS-1239` = "LincsRnaSeqSigData",
		      LINCSRS = "LincsRnaSeqSigData",
		      GDS = "GdsSigData",
		      CMAP = "CmapSigData",
		      ENC = "EncodeSigData",
		      LINCSTP = "LincsP100SigData",
		      LINCSOE = "LincsOeSigData",
		      DM = "DrugMatrixSigData",
		      CTRS = "ctrsSigData",
		      EBI = "ebiSigData",
          PG = "pharmGenSigData"
  )
# 	pval <- ifelse(sigType=="ENC", "prob", "P_Value")

	tmp <- DBI::dbGetQuery(mycon, paste0("select count(distinct geneID) from ",sigTable," where signatureID='",sigID,"' and geneID in ('", glist, "');"))[1,1]
	return(tmp)	
}
