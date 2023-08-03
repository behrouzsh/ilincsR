
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
#' @keywords test
#' @export 
#' @examples
#' res <- batchTregGRS(....)

filter <- function(filterchk,filters,parsedSData,includeORexclude){
tmpPSD<-NULL
if( !is.null(filterchk) && filterchk!="n"){
	fs<-unlist(strsplit(filters,","))
	uniqueProps<-NULL
	for(f in fs){
		pt<-unlist(strsplit(f,":"))[1]
		vl<-unlist(strsplit(f,":"))[2]
		if(includeORexclude == "1")
		{
			if(is.null(tmpPSD))
			{
				tmpPSD<-rownames(parsedSData[parsedSData[,pt]==vl,])
				uniqueProps[1]<-pt
			}
			else
			{
				if(is.element(pt,uniqueProps))
				{
					tmpPSD<-union(tmpPSD,rownames(parsedSData[parsedSData[,pt]==vl,]))
				}
				else
				{
					tmpPSD<-intersect(tmpPSD,rownames(parsedSData[parsedSData[,pt]==vl,]))
					uniqueProps[length(uniqueProps) +1]<-pt
				}
			}

		}
		else
		{
			parsedSData<-parsedSData[parsedSData[,pt]!=vl,]
		}

	}
	if(includeORexclude == "1")
	{
		parsedSData<-parsedSData[tmpPSD,]
	}
}
return(parsedSData)


}
