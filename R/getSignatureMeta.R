
#' A function to extract signature meta data from ilincs.org portal
#'
#' This function allows you to download/extract selected meta data related to specified signatures from ilincs portal.
#' @param prop A string containing a single or list of meta data columns from following list separated by commas: 
#'	c(libraryID, signatureID, datasetID, platform, factor, Level1, Level2, cID, compound, concentration, cellLine, 
#'	n_trt_samples, n_ctr_samples, antibodyTarget, peakType, treatment, time, lincsPertID, concordanceTable, stitchID,
#'	perturbagenID, pubChemID, treeviewSigId)
#'	By default the function returns 11 most important property of signatures.
#' @param signatures A string containing a single or list of signatures exist in ilincs portal separated by commas.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param author M. Fazel-Najafabadi
#' @keywords Signature.
#' @export 
#' @examples
#' ## not run
#' mysig_property <- getSignatureMeta(signatures="LINCSSH_1,LINCSSH_10,LINCSSH_100,LINCSSH_10000,ENC_112,GDS_152")
#' ## end not run
getSignatureMeta <- function(signatures, prop="signatureID,datasetID,perturbagenID,compound,lincsPertID,integratedGeneTargets,treatment,concentration,cellLine,time,factor,Level1,Level2,antibodyTarget", 
				    debuging=FALSE, sigDB="ilincs_sigs") 
{
	signatures <- gsub(",", "','", signatures)
	if(is.null(prop)) prop <- "signatureID"
	if(is.na(prop)) prop <- "signatureID"
	if(prop=="") prop <- "signatureID"
	servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
        mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=servSet$sigDB, host=servSet$host, port = servSet$port, password="public")
        sql <- paste0("select signatureID,", prop, " from signaturesMeta3 where signatureID in ('", signatures, "')")
        meta <- DBI::dbGetQuery(mycon, sql)
        DBI::dbDisconnect(mycon)
        colnames(meta) <- c("signatureID", unlist(strsplit(prop, ",")))
        return(meta)
}
