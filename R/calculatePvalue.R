
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

calculatePvalue <- function(pData,ExprData,prop,onetail)
{
	a<-Sys.time()
	ord <- order(pData[,prop])
	pData <- pData[ord,]
	#NoOfProbes <- dim(Values)[1]/dim(pData)[1]
	#geneNames <- Values[1:NoOfProbes,"Gene"]
	
#	ExprData <- NULL
#	ExprData<-matrix(Values[,"Value"],ncol=dim(pData)[1])	
#	unqVals <- as.character(unique(Values[,1]))
#	unqIndex <- match(unqVals,rownames(pData))
	
#	rownames(ExprData) <- as.vector(geneNames)
#	ExprData<-as.data.frame(ExprData)

	Pvalues<-NULL	
	## check missing data
        countingMissingData<-apply(ExprData,1,function(x) sum(is.na(x)))
        nonMissingData<-(countingMissingData != length(colnames(ExprData)))
        ExprData<-ExprData[nonMissingData,]
        ExprData<-as.matrix(ExprData)
#	rownames(ExprData) <- as.vector(geneNames[nonMissingData])

	status <-as.factor(as.character(pData[,prop]))
        design <- model.matrix(~ status)
	fit <- limma::lmFit(ExprData, design = design)
	fit2<-NULL
   	try(fit2 <- limma::eBayes(fit),TRUE)

	if(is.null(fit2))	
	{
		return(paste("Error in anova",geterrmessage()))
	}
	
	## if exactly 2 levels use p.values else use F.p.value
	if(length(unique(pData[,prop])) == 2) {	
		Pvalues <- fit2$p.value[,2]
	}
	else {
		Pvalues <- fit2$F.p.value
		names(Pvalues)<-names(fit2$p.value[,2])
	}
 	upPvalues <-NULL
	dnPvalues <- NULL 
	if(onetail)
	{
		upPvalues<- Pvalues/2
		names(upPvalues)<-names(Pvalues)
		dnPvalues<- 1 - Pvalues/2
		names(dnPvalues)<-names(Pvalues)
	}

	retlist<-NULL
	retlist<-list(Pvalues,upPvalues,dnPvalues)

	b<-Sys.time()
	print(paste("Start time for pval",a))
	print(paste("End time pval",b))
	print(paste("time for pValue",b-a))


	return(retlist)	
}
