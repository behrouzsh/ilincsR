
#' A function for getting ilincs portal information.
#'
#' This function allows you to see what R version, ilincsR version, compute server and user is being used in the analysis.
#' @param author M. Fazel-Najafabadi
#' @keywords ...
#' @export
#' @examples
#' ilincsInfo()

ilincsInfo <- function() {
#     require(ilincsR)
    ses <-  sessionInfo()
    sys <- as.list(Sys.info())
    var <- Sys.getenv()
    class(var) <- "vector"
    var <- var[names(var) %in% c("chost", "dhost", "dport", "sigDB")]
    for(v in c("chost", "dhost", "dport", "sigDB")) if (is.na(var[v])) var[v] <- NA

    res <- list(
		R.version = ses$R.version$version.string,
		ilincsRpkg = paste(ses$otherPkgs[[grep("ilincsR", ses$otherPkgs)]]$Package, ses$otherPkgs[[grep("ilincsR", 
				    ses$otherPkgs)]]$Version, sep="_"),
		computeServer = sys$nodename,
		dbHost = var[["dhost"]],
		dbPort = var[["dport"]],
		sigDB = var[["sigDB"]],
		user = sys$user
		)
    res
}
