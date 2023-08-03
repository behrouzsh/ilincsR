
#' A function to extract gene targets for a drug.
#'
#' This function allows you to extract gene targets related to specified drugID, lsmID or signatureID from ilincs portal.
#' @param lsmID  A vector of LSM IDs. 
#' @param drugID A vector of BRD IDs (perturbagenIDs).
#' @param signatureID A vector of chemical perturbage (CP) signature IDs.
#' @param sigDB  Signature database to use. 
#' @param debuging For debuging purposes only 
#' @param author M. Fazel-Najafabadi
#' @keywords Signature, perturbage.
#' @export 
#' @examples
#' ## not run
#' targets <- drugTargets(lsmID=c("LSM-45711", "LSM-5242", "LSM-45954"))
#' ## end not run

drugTargets <- function(lsmID, drugID=NULL, signatureID=NULL, sigDB="ilincs_sigs", debuging=FALSE) {
  servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=servSet$sigDB, host=servSet$host, port = servSet$port, password="public")
  
  if(!is.null(signatureID)){
    sigs <- paste(unique(signatureID), collapse = "','")
    # sql <- paste0("select distinct meta.signatureID,tar.pertID,tar.lsmID,tar.geneSymbol from  LincsCpGeneTargets tar, signaturesMeta3 meta where meta.libraryID='LIB_5' and meta.signatureID in('", sigs, "') and pert.pertID in meta.perturbagenID")
    sql <- paste0("select distinct meta.signatureID, meta.perturbagenID from signaturesMeta3 meta where meta.libraryID='LIB_5' and meta.signatureID in('", sigs, "')")
    drugSig <- DBI::dbGetQuery(mycon, sql)
    rownames(drugSig) <- drugSig$signatureID
    drugSig <- drugSig[signatureID,]
    drugs <- paste(unique(drugSig$perturbagenID), collapse = "','")
    sql <- paste0("select distinct pertID,lsmID,geneSymbol from  LincsCpGeneTargets where pertID in ('", drugs, "')")
    targets <- DBI::dbGetQuery(mycon, sql)
    
    spl <- split(targets, targets$pertID)
    spl <- sapply(spl, function(j) paste(j$geneSymbol, collapse = ", "))
    mch <- match(drugSig$perturbagenID, names(spl))
    res <- cbind.data.frame(drugSig, targets=spl[mch], stringsAsFactors = F)
    res <- res[,-2]
  } else if(!is.null(drugID)) {
    drugs <- paste(unique(drugID), collapse = "','")
    sql <- paste0("select distinct pertID,lsmID,geneSymbol from  LincsCpGeneTargets where pertID in ('", drugs, "')")
    targets <- DBI::dbGetQuery(mycon, sql)
    spl <- split(targets, targets$pertID)
    spl <- sapply(spl, function(j) paste(j$geneSymbol, collapse = ", "))
    mch <- match(drugID, names(spl))
    res <- data.frame(drugID = drugID, targets=spl[mch], stringsAsFactors = F)
    
  } else {
    lsms <- paste(unique(lsmID), collapse = "','")
    sql <- paste0("select distinct pertID,lsmID,geneSymbol from  LincsCpGeneTargets where lsmID in ('", lsms, "')")
    targets <- DBI::dbGetQuery(mycon, sql)
    spl <- split(targets, targets$lsmID)
    spl <- sapply(spl, function(j) paste(j$geneSymbol, collapse = ", "))
    mch <- match(lsmID, names(spl))
    res <- data.frame(lsmID = lsmID, targets=spl[mch], stringsAsFactors = F)
  }
  DBI::dbDisconnect(mycon)
  return(res)
}
