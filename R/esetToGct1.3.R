
#' A function to convert an eset (expression set class) to a "gct" file format
#'
#' This function allows you to convert an eset (Biobase class object) to a "gct" format (text file) and save the ".gct" file.
#' @param eset The expression set (expression set class from Biobase).
#' @param gctFile This optional parameter specifies the resulting "gct" file name. If not provided, the experiment name will be used.
#' @param path_to_write This parameter specifies the path to save the "gct" format of the downloaded data.
#' @param author M. Fazel-Najafabadi
#' @keywords gct
#' @export 
#' @examples
#' ## not run
#' mygct <- esetToGct1.3(eset=eset, gctFile="my_gct_file", path_to_write="your/path/to/write/here")
#' ## end not run

esetToGct1.3 <- function(eset, gctFile=NA, path_to_write) {
#	require(Biobase)
#	assign(esetName, get(load(eset)))
	edata <- Biobase::exprs(eset)
	fdata <- Biobase::fData(eset)
	pdata <- Biobase::pData(eset)
	endChar <- function (x, suffix) {
	    if (!is.character(x) || !is.character(suffix)) stop("non-character object(s)")
	    n = nchar(x)
	    suppressWarnings(substr(x, n - nchar(suffix) + 1L, n) == suffix)
	}
	
	fdataFrame <- as.data.frame(matrix("na", nrow=dim(pdata)[2], ncol=dim(fdata)[2]), stringsAsFactors=FALSE)
	names(fdataFrame) <- names(fdata)
	fdataFrame <- rbind(fdataFrame, fdata)
	df <- cbind(fdataFrame, rbind(t(pdata), edata), stringsAsFactors=FALSE)
	df <- rbind(colnames(df), df)
	df <- cbind(c("id", colnames(pdata), rownames(fdata)), df, stringsAsFactors=FALSE)
	colnames(df) <- NULL
	rownames(df) <- NULL

	path_to_write <- sub("/?$", "", path_to_write)
	path_to_write <- sub("$", "/", path_to_write)
	if (is.na(gctFile)) {
#		gctFile <- paste(gctDir, sub("_*eset\\.RData$", "\\.gct", deparse(substitute(esetName))), sep="")
		print(deparse(substitute(eset)))
		gctFile <- paste(path_to_write, deparse(substitute(eset)), ".gct", sep="")
#		print(gctFile)
	} else {
		gctFile <- ifelse(endChar(gctFile, ".gct") , paste0(path_to_write, gctFile), paste0(path_to_write, gctFile, ".gct"))
	}
	write(paste("#1.3\n", dim(edata)[1], "\t", dim(edata)[2], "\t", dim(fdata)[2], "\t", dim(pdata)[2], sep=""), file=gctFile)
	write.table(df, file=gctFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}
