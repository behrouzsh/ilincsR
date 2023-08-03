
#' A function for generate session identifier.
#'
#' This function allows you to generate a unique session ID at random.
#' @param author M. Fazel-Najafabadi
#' @keywords ...
#' @export
#' @examples
#' generateSID()

generateSID <- function() {
    
    sessionID <- paste(c(strsplit(date(),split=" ")[[1]],as.integer(runif(1)*10e6)),sep="",collapse="_")
    sessionID <- gsub(":","_", sessionID)
    sessionID <- gsub("__","_", sessionID)
    return(sessionID)

}
