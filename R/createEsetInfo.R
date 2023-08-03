
#' A function for ...
#'
#' This function allows you to ...
#' @param exp This is basically the experiment name.
#' @param path_to_write This parameter specifies the path which user wants to save the results.
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @export
#' @examples
#' res <- createEsetInfo(....)

createEsetInfo <- function(exp, glist, path_to_write="/mnt/raid/tmp/", up=4000, down=1000, window.size=50, debuging=FALSE) {

# genes1 <- find_genes_hgenes_in_platform(exp, glist, test)
# 
# if (is.null(genes1) || genes1=="") {
#     # return(list(c("No genes found", "NA", "NA"), NULL))
#     message("Non of the genes found in the platform")
#     return(NULL)
#   }
pre <- preprocessing(glist, exp, path_to_write, up, down, window.size, debuging=debuging)
if (pre[[1]] == "No genes found") {
    message("Non of the genes found in the platform")
    return(list(esetName = "Non of the genes found in the platform", 
	    foundGenes = NA, 
	    NoOfGenes = NA, 
	    NoOfProbes = NA, 
	    NoOfSamples = NA, 
	    xls = NA))
}
genes <- unlist(pre)[3]
esetName <- unlist(pre)[1]

#genes <- gsub("NA,","",genes)
#genes <- gsub(",NA","",genes)
#genes <- gsub(",,","",genes)

genes <- gsub("NA", "", genes)
glist <- na.omit(as.integer(unlist(strsplit(genes, ','))))
NoOfGenes <- length(glist)

esetPath <- paste0(path_to_write, esetName, ".RData")
tmp <- load(esetPath)
eset <- get(tmp)
dims <- dim(Biobase::exprs(eset))
  #noProbes<-dim(Biobase::exprs(esetName))[1]

#com <- paste0("dim(Biobase::exprs(", esetName, "))[1]")
#NoOfProbes <- eval(parse(text=com))
NoOfProbes <- dims[1]

#com <- paste0("dim(Biobase::exprs(", esetName, "))[2]")
#NoOfSamples <- eval(parse(text=com))
NoOfSamples <- dims[2]

#com <- paste0("create_xls_from_eset(", esetName, ",path_to_write)")
#xlsID <- eval(parse(text=com))
#xls <- paste0("http://www.ilincs.org/tmp/download_data", xlsID, ".xls")
xlsID <- create_xls_from_eset(eset, path_to_write)
## esetToGct1.3(eset, paste0("/filteredeset_", sessionID), path_to_write)
if (path_to_write == "/mnt/raid/tmp/") {
#     xls <- paste0("http://www.ilincs.org/tmp/download_data", xlsID, ".xls")} else {xls <- paste0(path_to_write, xlsID, ".xls")
    xls <- paste0("download_data", xlsID, ".xls")} else {xls <- paste0(path_to_write, xlsID, ".xls")
}

res <- list(esetName = esetName, 
	    foundGenes = glist, 
	    NoOfGenes = NoOfGenes, 
	    NoOfProbes = NoOfProbes, 
	    NoOfSamples = NoOfSamples, 
	    xls = xls)
return(res)
}
