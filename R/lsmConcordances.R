
#' A function to retrieve drug similarity
#'
#' This function allows you to retrieve necessary metadata about an experiment from databases.
#' @param lsmList A vector of LSM IDs for drugs.
#' @param minSimilarity A single numeric value for minimum similarity between drugs, default to "0".
#' @param author M. Fazel-Najafabadi
#' @keywords LSM ID.
#' @export 
#' @examples
#' ## not run
#' myLSMs <- lsmConcordances(c("LSM-10001","LSM-10000","LSM-36603","LSM-1535"), minSimilarity=0.2, debuging=T)
#' ## end not run

lsmConcordances <- function(lsmList, minSimilarity=0, debuging=FALSE) {
    servSet <- getServerSettings(debuging)
    if(length(lsmList) > 1) lsmList <- paste(lsmList, collapse="','")
    mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "ChemicalDB", port = servSet$port, host = servSet$host, password = "public")
    lsms <- DBI::dbGetQuery(mycon, paste0("select firstChemID,secondChemID,similarity from TanimotoIndex where firstChemID in ('", lsmList, 
					  "') and secondChemID in ('", lsmList, "') and similarity >= ", minSimilarity))
    DBI::dbDisconnect(mycon)
    return(lsms)
}

