
#' A function for finding common genes between a list and an experiment.
#'
#' This function allows you to know how many of the genes in your list are or have homologs in a spesified experiment.
#' @param exp ...
#' @param glist ...
#' @param test ...
#' @param authors M. Fazel-Najafabadi
#' @keywords test
#' @export 
#' @examples
#' ## Not run
#' genelist <- "1002,9447,10001,10096,79575"
#' res <- myGenesInExp(exp="EDS-1014", glist=genelist, debuging=FALSE)
#' ## End, not run

myGenesInExp <- function(exp, glist, org="Hs", debuging=FALSE) {
	tmp <- find_genes_hgenes_in_platform(exp, glist, org=org, debuging=debuging)
	tmp <- length(unlist(strsplit(tmp, ",")))
	return(tmp)	
}
