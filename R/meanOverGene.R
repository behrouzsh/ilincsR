
#' A function to average expression and p-Values over unique genes.
#'
#' This function allows you to create a signature with unique gene IDs by averaging over multiplicated ones.
#' @param rawSignature A data frame of three or four columns (gene ID, gene Name, coefficients and possibly p-Values).
#' @param author M. Fazel-Najafabadi
#' @keywords ...
#' @examples
#'	## not run
#'	mysig <-  meanOverGene(mysig)
#'	## end not run

meanOverGene <- function(rawSignature) {
	
	if (dim(rawSignature)[2] == 3) {
	    colnames(rawSignature) <- c("geneID", "geneName", "coefficients")
	} else {
	    colnames(rawSignature) <- c("geneID", "geneName", "coefficients", "Pvals")
	}
	geneID <- as.vector(rawSignature$geneID)
	geneName <- as.vector(rawSignature$geneName)
	Signature <- aggregate(rawSignature[,-(1:2),drop=FALSE], by = list(geneID, geneName), FUN = "mean", na.rm=TRUE)
	colnames(Signature)[1:2] <- c("geneID", "geneName")
# 	if(dim(rawSignature)[2] == 4) Signature <- aggregate(rawSignature[,-(1:2),drop=FALSE], by = list(geneID, geneName), FUN = "mean", na.rm=TRUE)
	
# 	if(dim(rawSignature)[2] == 4) unqIdPvals <- split(rawSignature[,"Pvals"], rawSignature[,"geneID"])
# 	if(dim(rawSignature)[2] == 3) Signature <- data.frame(geneID=names(unqIdCoeffs), coefficients=mergeCoeffs, stringsAsFactors=FALSE)
# 	if(dim(rawSignature)[2] == 4) Signature <- data.frame(geneID=names(unqIdCoeffs), coefficients=mergeCoeffs, Pvals=mergePvals, stringsAsFactors=FALSE)
	return(Signature)
}
