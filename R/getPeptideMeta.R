
#' A function to retrieve peptide info
#'
#' This function allows you to retrieve peptide metadata from uniprot.
#' @param peptideList A vector of LSM IDs for drugs.
#' @param metacols Metadata column names needed.
#' @param debuging Only for database test purposes.
#' @param author M. Fazel-Najafabadi
#' @keywords Uniprot.
#' @export 
#' @examples
#' ## not run
#' myLSMs <- getPeptideMeta(c("LSM-10001","LSM-10000","LSM-36603","LSM-1535"), minSimilarity=0.2, debuging=T)
#' ## end not run

getPeptideMeta <- function(peptideList, metacols="*", debuging=FALSE) {
    servSet <- getServerSettings(debuging)
    if(length(peptideList) > 1) peptideList <- paste(peptideList, collapse="','")
    mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "ProteinDB", port = servSet$port, host = servSet$host, password = "public")
    meta <- DBI::dbGetQuery(mycon, paste0("select ", metacols, " from uniprot where accession in ('", peptideList, "')"))
    DBI::dbDisconnect(mycon)
    return(meta)
}

