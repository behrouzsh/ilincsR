
#' A function to convert Gene IDs to Gene Symbols
#'
#' This function allows you to convert official Gene IDs (Entrez) to Gene Symbols.
#' @param idlist A character string containing a single or list of Entrez gene IDs separated by commas.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param species Three different available species in the Genomics Portal database ("Hs" for Human, "Mm" for Mouse or "Rn" for Rat).
#' @keywords Entrez Gene ID
#' @export 
#' @examples
#' geneSymbols <- geneid2symbol(idlist = "57535,11010,283987,2920,2113,10612,9619")

geneid2symbol <- function(idlist, species="Hs", description=FALSE, verbose=FALSE, debuging=FALSE) {

##tax <- org2tax(species)
servSet <- getServerSettings(debuging=debuging)
mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDB", port = servSet$port, host = servSet$host, password = "public")
#if (!test) {
#	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", port = 4040,host="10.165.4.231", password='public')
#	} else {
#rs <- DBI::dbSendQuery(mycon, sql)
#gene_table <- DBI::fetch(rs, n = -1)

#@sql <- paste0("select GeneID,Symbol from GeneInfo where GeneID=", idlist[id], " and TaxID=", tax)
idlist <- gsub(",", "','", idlist)
ids <- (unlist(strsplit(idlist, "','")))#unique
sql <- paste0("select GeneID,Symbol from GeneInfo where GeneID in ('", idlist, "')")
if(description) sql <- paste0("select GeneID,Symbol,Description from GeneInfo where GeneID in ('", idlist, "')")
#print(paste(sql))
gene_table <- data.frame(GeneID=ids, Symbol=NA, stringsAsFactors=FALSE)# NULL
tmp1 <- DBI::dbGetQuery(mycon, sql)
mm <- match(gene_table$GeneID, tmp1$GeneID)

notfounds <- ids[(ids %in% tmp1$GeneID)==FALSE]
    if (length(notfounds) > 0 & verbose) {
	print(paste0("Symbol for gene ID(s) '", paste(notfounds, collapse=", "),"' not found"))
    }
gene_table[, "Symbol"] <- tmp1[mm, "Symbol"]
if(description) gene_table$Description <- tmp1[mm, "Description"]

DBI::dbDisconnect(mycon)
return(na.omit(gene_table))
}
