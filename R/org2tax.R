
#' A function for converting organism name to tax ID.
#'
#' This function allows you to convert organism character code to tax ID.
#' @param org Character value for organism.
#' @param authors M. Fazel-Najafabadi
#' @keywords test
#' @examples
#' ## Not run
#' tax <- org2tax("Hs")
#' ## End, not run

org2tax <- function(org) {
    if(org %in% c("Hs", "human")) tax <- 9606 else if(org %in% c("Mm", "mouse")) tax <- 10090 else if(org %in% c("Rn", "rat")) tax <- 10116 else stop("Not appropriate organism!")
#     tax <- switch(org, Hs = 9606, Mm = 10090, Rn = 10116)
    return(tax)	
}
