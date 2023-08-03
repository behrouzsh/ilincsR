
#' A function for choosing database and coresponding server.
#'
#' This function allows you to choose a server which provides databases. It is for testing purposes only before a dataset/database is released.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param sigSet Which signature version user wants to use in the analysis.
#' @param author M. Fazel-Najafabadi
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

getServerSettings <- function(debuging=FALSE, sigDB="ilincs_sigs", chost="dev.ilincs.org", dhost="gimm2.ketl.uc.edu", dport=3306) 
{
# 	if (test == TRUE & sigSet=="new")  {return(data.frame(sigDB = "ilincs_new", host = "gimm2.ketl.uc.edu", port = 3306, stringsAsFactors = FALSE))}
# 	if (test == TRUE & sigSet=="old")  {return(data.frame(sigDB = "ilincsTest2", host = "gimm2.ketl.uc.edu", port = 3306, stringsAsFactors = FALSE))}
# 	if (test != TRUE & sigSet=="new") {return(data.frame(sigDB = "ilincs_new", host = "eh3.uc.edu", port = 4040, stringsAsFactors = FALSE))} 
# 	if (test != TRUE & sigSet=="old") {return(data.frame(sigDB = "ilincsTest2", host = "eh3.uc.edu", port = 4040, stringsAsFactors = FALSE))}
	
	if(debuging) {
	    host <- dhost
	    compute <- chost
	    port <- dport
	    sigDB <- sigDB
	} else {
	    vars <- Sys.getenv()
	    class(vars) <- "vector"
	    vars <- vars[names(vars) %in% c("chost", "dhost", "dport", "sigDB")]
	    for(v in c("chost", "dhost", "dport")) if (is.na(vars[v])) return(print(paste0("Missing ", v)))
	    host <- vars["dhost"]; compute <- vars["chost"]; port <- as.integer(vars["dport"]); sigDB <- vars["sigDB"]
	}
	
	if (sigDB=="ilincs_sigs")  {return(data.frame(sigDB = "ilincs_sigs", compute=compute, host = host, port = port, stringsAsFactors = FALSE, row.names=NULL))}
	if (sigDB=="ilincs_new")  {return(data.frame(sigDB = "ilincs_new", compute=compute, host = host, port = port, stringsAsFactors = FALSE, row.names=NULL))}
	if (sigDB=="ilincsTest2") {return(data.frame(sigDB = "ilincsTest2", compute=compute, host = host, port = port, stringsAsFactors = FALSE, row.names=NULL))}
	
	
#   # return(data.frame( sigDB='ilincsTest2', host='db3.ketl.uc.edu', port = '3306'))
#   # return(data.frame( sigDB='ilincsTest2', host='gimm5.ketl.uc.edu', port = '3306'))
}
