
#' A function to download a nicely annotated whole dataset from ilincs.org portal
#'
#' This function allows you to download a dataset as "gct" format with gene annotations.
#' @param exp The experiment which you want to filter and download. It should be present in the portal's database.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param author M. Fazel-Najafabadi
#' @keywords Download.
#' @export 
#' @examples
#' ## not run
#' mygct <- getAnnotatedExp(exp="EDS-1014")
#' ## end not run

getAnnotatedExp <- function(exp, id.na.rm=TRUE, debuging=FALSE) {
    servSet <- getServerSettings(debuging)
    mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
    esetinfo <- DBI::dbGetQuery(mycon, paste0("select Platform, eset_path, Organism from ExperimentMetadata where experiment = '", exp, "'"))
    species <- esetinfo$Organism
    tax <- switch(species, human = "Hs", Human = "Hs", mouse = "Mm", Mm = "Mm", rat = "Rn")
    platform <- esetinfo$Platform
    platformID <- DBI::dbGetQuery(mycon, paste0("select ID from Platforms where Platform = '", platform, "'"))[1,1]
    esetFileName <- esetinfo$eset_path
    esetFileName <- gsub("/mnt/", "/opt/", esetFileName)
    eset <- get(load(esetFileName))
    pg <- DBI::dbGetQuery(mycon, paste0("select Probe, Gene from ProbeGene where Platform_ID = '", platformID, "'"))
    DBI::dbDisconnect(mycon)
    probe <- rownames(Biobase::exprs(eset))
    geneID <- pg$Gene[match(probe, pg$Probe)]
    idName <- geneid2symbol(paste(geneID, collapse=","), species=tax, debuging=debuging, description=TRUE)
    mm <- match(geneID, idName$GeneID)
    fdata <- data.frame(PROBE=probe, Name_GeneSymbol=idName$Symbol[mm], ID_geneid=as.character(idName$GeneID[mm]), DESCRIPTION=idName$Description[mm], stringsAsFactors=FALSE)
    rownames(fdata) <- fdata$PROBE
    if(id.na.rm) fdata <- fdata[!is.na(fdata$ID_geneid), ]
    eset <- eset[rownames(fdata),]
    Biobase::fData(eset) <- fdata
    tp <- apply(Biobase::pData(eset),2,class)
    tq <- which(tp=="list")
    if (length(tq)>0) {
	pdata <- Biobase::pData(eset)
	for(c in tq) pdata[,c] <- as.character(pdata[,c])
	Biobase::pData(eset) <- pdata
    }
    return(eset)
}


