
#' A function to convert Gene Symbols to Gene IDs
#'
#' This function allows you to convert official Gene Symbols to Gene IDs (Entrez gene IDs).
#' @param symlist A character string containing a single or list of gene symbols separated by commas.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param exact Logical, to find exactly matched gene symbols to database. If FALSE, matching is not case sensitive.
#' @param species Three different available species in the Genomics Portal database ("Hs" for Human, "Mm" for Mouse or "Rn" for Rat).
#' @keywords Entrez Gene ID
#' @export 
#' @examples
#' geneSymbols <- geneid2symbol(idlist = "57535,11010,283987,2920,2113,10612,9619")


symbol2geneid <- function(symlist, org="Hs", debuging=FALSE) {

tax <- org2tax(org)
symlist <- gsub(",", "','", symlist)
servSet <- getServerSettings(debuging)
mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDB", port = servSet$port, host = servSet$host, password = "public")
#if (!test) {
#	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", port = 4040,host="10.165.4.231", password='public')
#	} else {
#	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", host="gimm2.ketl.uc.edu", password='public')}
#rs <- DBI::dbSendQuery(mycon, sql)
#gene_table <- DBI::fetch(rs, n = -1)

# syms <- strsplit(symlist,",")
# syms <- paste(syms, collapse="','")
syms <- unlist(strsplit(symlist, "','"))
sql <- paste0("select GeneID,Symbol from GeneInfo where Symbol in ('", symlist, "') and TaxID=", tax)
## # # # # sql <- paste0("select GeneID,Symbol,TaxID from GeneInfo where Symbol in ('", symlist, "')")
tmp <- DBI::dbGetQuery(mycon, sql)
#@!}
gene_table <- tmp
#@!!!gene_table <- gene_table[!duplicated(gene_table$Symbol), , drop=FALSE]
#DBI::dbClearResult(RMySQL::dbListResults(mycon)[[1]])
DBI::dbDisconnect(mycon)
notFounds <- syms[!(syms %in% tmp$Symbol)]
if (length(notFounds) > 0) print(paste("Entrez ID for gene symbols '", paste(notFounds, collapse=","), "' not found"))
#@!if (exact) gene_table <- gene_table[(gene_table$Symbol %in% strsplit(symlist, ",")), ] ## to find exact match symbols to database
return(gene_table)
}


