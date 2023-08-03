
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

static_heatmaps <- function(sessionID,dataClust,nFactorLevels,factorLevels,geneNames,pal,rsd,allColors,prop,nColors,
			    ifCluster,prop_distance,isPropNone,path_to_write,dataFormat, sigDB="ilincs_sigs", chost="ilincs.org", debuging=FALSE)
{
##    library(stringr)
	if(dim(dataClust)[1]==1) {
		dataClust <- rbind(dataClust[1,],dataClust[1,])
	}

	rowLable <- apply(as.matrix(dataClust[,2]),1, function(x) paste(unlist(strsplit(x,":"))[1], sep = ":"))
        rowLableProbe <- apply(as.matrix(dataClust[,2]),1, function(x) paste(unlist(strsplit(x,":"))[2], sep = ":"))
        ## Mehdi
        #print(paste("**********",  max(nchar(rowLableProbe))))
        #rowLableProbe <- as.vector(sapply(rowLableProbe, function(x) paste(x, paste(rep(" ", 15+max(nchar(rowLableProbe))-nchar(x)), 
	#							collapse="", sep=""), collapse="", sep=""))) ## I wraped the names 
	#as.vector(sapply(string, function(x) paste(x, paste(rep(" ", max(nchar(string))-nchar(x)), collapse=""), sep="")))
        #rowLable <- paste(rowLable, rowLableProbe, sep = ":")
        #rowLableProbe <- str_pad(rowLableProbe, width=max(nchar(rowLableProbe))+1, pad=" ", side="right") ## I wraped the names 
        #rowLable <- paste(rowLable, rowLableProbe, sep = " : ")
        rowLable <- paste(rowLableProbe, rowLable, sep = " : ") ## I switched the gene and ensemble Id
	perGene <- 180/1000
        heightRatio <- 200/1000
        totalCols <- length(colnames(dataClust))
	currRowOrder <- order(rowLable)
	##Original data
	pdf(paste(path_to_write,"temp",sessionID,".pdf",sep=""),height=max(13,(perGene*length(geneNames))), width= 11)

	colnam <- colnames(dataClust)
	if(dim(dataClust)[1]==1) {
		dataClust <- rbind(dataClust[1,],dataClust[1,])
		rowLable <- rbind(rowLable,"")
		currRowOrder <- order(rowLable)
	}
	colnames(dataClust) <- colnam
	colsepindex <- as.vector(seq(1,dim(dataClust[,-(1:2)])[2],1))
	rowsepindex <- as.vector(seq(1,dim(dataClust[,-(1:2)])[1],1))
	
	cluster_datasets(sessionID,dataClust,nFactorLevels,ifCluster,currRowOrder,rowLable,rowsepindex,colsepindex,rsd,pal,prop_distance,isPropNone,FALSE,
			path_to_write, prop, factorLevels, allColors, nColors, sigDB=sigDB, chost=chost, debuging=debuging)
	## Mehdi: I added one line to test!
	#legend(x=0.1,y=1,title = prop,legend=factorLevels,col=allColors[(((1:nFactorLevels)-1)%%nColors)+1],pch=rep(19,nFactorLevels),ncol=3,bty="o",bg = "lightgrey",pt.cex=2)

	try(dev.off()) ## Mehdi: need to move this part to "cluster_datasets"
#@@	if(!isPropNone)
#@@	{
#@@	legendHeight = (nFactorLevels/3)*0.3
#@@	pdf(paste(path_to_write,"legendtemp",sessionID,".pdf",sep=""), height = max(8.5,legendHeight), width = 11)
#@@	plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",cex.main=2)
#@@	legend(x=0.1,y=1,title = prop,legend=factorLevels,col=allColors[(((1:nFactorLevels)-1)%%nColors)+1],pch=rep(19,nFactorLevels),ncol=3,bty="o",bg = "lightgrey",pt.cex=2)
#@@	try(dev.off())
#@@	}

	##Centered data
	if(dataFormat != "Gff") {
	means <- apply(data.matrix(dataClust[,-(1:2)]),1,median,na.rm=T)
	centeredDataClust <- sweep(data.matrix(dataClust[,-(1:2)]),1,means,"-")
	centeredDataClust <- cbind(dataClust[,1:2],centeredDataClust)

#	source("http://eh3.uc.edu/r/gimmHeat.R")
	pdf(paste(path_to_write,"tempgheat",sessionID,".pdf",sep=""),height=max(10,(perGene*length(geneNames))), width=11)

	cluster_datasets(sessionID,centeredDataClust,nFactorLevels,ifCluster,currRowOrder,rowLable,rowsepindex,colsepindex,rsd,pal,prop_distance,isPropNone,TRUE,
			path_to_write, prop, factorLevels, allColors, nColors, sigDB=sigDB, chost=chost, debuging=debuging)

	try(dev.off())
	}
}
