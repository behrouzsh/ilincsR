
#' A function for ...
#'
#' This function allows you to ...
#' @param query.p ...
#' @param reference.p ...
#' @param query.gL ...
#' @param reference.gL ...
#' @param na.rm ...
#' @param estimateNullDistr ...
#' @param nullDistrQuantiles ...
#' @param nullDistrN ...
#' @param tolerateWarnings ...
#' @param pAdjust.method.query ...
#' @param pAdjust.method.reference ...
#' @param lambda ...
#' @param rescale ...
#' @param plotRescaled ...
#' @param multicoreN ...
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

get_genepairs_freq <- function(gp){
	rows <- strsplit(gp,":")[[1]]
        cols <- NULL;datasets <- NULL;table<-NULL;
        for(a in 1:length(rows)){
                cols <- unlist(strsplit(rows[a],","))
		if(cols[4] == "null") next
                g1<-cols[1]
                g2<-cols[2]
                freq<-cols[3]
                for(b in 4:length(cols)){
                        if(is.null(datasets)){
                                 datasets <- paste(cols[b],sep=",")
                        } else {
                                datasets <- paste(datasets,cols[b],sep=",")
                        }
                }
                table <- rbind(table,data.frame(g1,g2,freq,datasets))
                datasets <- NULL
        }
}
