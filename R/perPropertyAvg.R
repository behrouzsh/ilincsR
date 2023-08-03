
#' A function for ...
#'
#' This function allows you to ...
#' @param dataClust ...
#' @param nFactorLevels ...
#' @param factorLevels ...
#' @param geneNames ...
#' @param pal,allColors,nColors ...
#' @param sessionID ...
#' @param dataType ...
#' @param path_to_write 
#' @param path_to_write This necessary parameter specifies the path which user wants to save the heatmaps output file.
#' @param Pvalues ...
#' @param PvalueCutoff ...
#' @param prop prop defines which property of the experiment is going to be used for statistical analysis to specify the groups of samples.
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @export
#' @examples ...
#' res <- perPropertyAvg()

perPropertyAvg <- function(dataClust, nFactorLevels, factorLevels, geneNames, pal, allColors, nColors, 
			   sessionID, dataType, path_to_write, Pvalues=NULL, PvalueCutoff=0.05, prop, debuging=FALSE, chost="ilincs.org") {
## Mehdi: for testing purposes
#save(dataClust, nFactorLevels, factorLevels, geneNames, pal, allColors, nColors, sessionID, dataType, 
#      path_to_write, Pvalues, PvalueCutoff, prop, file=paste0(path_to_write, "perProperty.RData"))
	perPropertyAverages <- NULL
	for(i in 1:nFactorLevels){
		temp <- which(as.vector(colnames(dataClust))==factorLevels[i])
		perPropertyAverages <- cbind(perPropertyAverages, apply(as.matrix(dataClust[,temp]),1,function(x) mean(x,na.rm=TRUE)))
	}
	colnames(perPropertyAverages) <- factorLevels
	fac <- as.factor(as.vector(colnames(dataClust)[-(1:2)]))

	if(is.null(Pvalues)){
		if(nFactorLevels == 1) #for one level sample
		{
			#t.test(ww)
			#anova(lm(ww~1))
			try(Pvalues <- apply(dataClust[,-(1:2)],1,function(x) t.test(x)$p.value))
		} else {
			try(Pvalues <- apply(dataClust[,-(1:2)],1,function(x) anova(lm(x~as.factor(colnames(dataClust)[-(1:2)])))[[5]][1]))
		}
		if(sum(is.na(Pvalues)) == length(Pvalues))
			return("error in anova")
	}
	retDf <- cbind(perPropertyAverages,Pvalues)
	
    servSet <- getServerSettings(debuging=debuging, chost=chost)
    if(servSet$compute == "dev.ilincs.org") {

	significanceBar <- rep("green",length(Pvalues))
	#significanceBar[Pvalues<0.05]<-"red"
	significanceBar[Pvalues<PvalueCutoff] <- "red"
	table(significanceBar)
	pvaluersd <- rep("blue",nFactorLevels)
	for(i in 1: nFactorLevels)
	{
		pvaluersd[i] <- allColors[((i-1)%%nColors)+1]
	}
	perPropertyAverages<-as.data.frame(perPropertyAverages)

	perGene <- 180/1000
	heightRatio <- 40/1000
	totalCols <- length(colnames(dataClust))
	
	rowLable <- apply(as.matrix(dataClust[,2]),1, function(x) paste(unlist(strsplit(x,":"))[1], sep = ":"))
	rowLableProbe <- apply(as.matrix(dataClust[,2]),1, function(x) paste(unlist(strsplit(x,":"))[2], sep = ":"))
	rowLable <- paste(rowLable,rowLableProbe, sep = ":")

	currRowOrder <- order(rowLable)
	try(perPropertyAverages <- as.data.frame(perPropertyAverages))

	pdf(paste(path_to_write,"pvaluetemp",sessionID,".pdf",sep=""),height=max(10,(perGene*length(geneNames))), width=11)

	if(dim(perPropertyAverages)[1]==1)
	{
		perPropertyAverages <- rbind(perPropertyAverages[1,],perPropertyAverages[1,])
		rowLable <- c(rowLable,"")
                currRowOrder <- rowLable
		significanceBar <- c(significanceBar,significanceBar[1])

	}
	if(dim(perPropertyAverages)[2]>1)
	{
		meanPerPropertyAverage <- apply(perPropertyAverages,1,mean,na.rm=T)
		for(i in 1:dim(perPropertyAverages)[2]) perPropertyAverages[,i]<-perPropertyAverages[,i]-meanPerPropertyAverage
	}

#	source("http://eh3.uc.edu/r/gimmHeat.R")
	rownames(perPropertyAverages) <- rowLable 
	colsepindex <- as.vector(seq(1,dim(perPropertyAverages)[2],1))

	if(dim(perPropertyAverages)[2]==1) #only 1 column
	{
		perPropertyAverages <- cbind(perPropertyAverages[,1],perPropertyAverages[,1]) # repeat column
		pvaluersd <- c(pvaluersd,pvaluersd[1])
		rownames(perPropertyAverages) <- rowLable
		colsepindex <- NULL
	}
	rowClust <- hclust(dist(data.matrix(perPropertyAverages)),method="complete")
        rowDendro <- as.dendrogram(rowClust)

	#try(gimmHeat(data.matrix(perPropertyAverages),gimmTrunc=c(-3,3),dendrogram="row",labRow=rownames(perPropertyAverages),labCol=NA,Colv=FALSE,Rowv=rowDendro,trace="none",scale = "none",col=pal,main="Statistical analysis of diffrence in average log2 levels",key = T,density="none",ColSideColors=pvaluersd, RowSideColors=significanceBar,margins = c(5,12),keysize = 0.8,lhei = c(1,max(2,(0.05 * length(geneNames)))),colsep=colsepindex,sepwidth=c(0.01,0.01),cexRow = min(0.3 + 1/log10(length(geneNames)),0.9)))
	try(gplots::heatmap.2(data.matrix(perPropertyAverages),dendrogram="row",labRow=rownames(perPropertyAverages),labCol=NA,Colv=FALSE,Rowv=rowDendro,trace="none",scale = "none",col=pal,main="Statistical analysis of difference in average log2 levels",key = T,density="none",ColSideColors=pvaluersd, RowSideColors=significanceBar,margins = c(8,25),keysize = 0.8,lhei = c(1,max(2,(0.05 * length(geneNames)))),colsep=colsepindex,sepwidth=c(0.01,0.01),cexRow = min(0.3 + 1/log10(length(geneNames)),0.9)))
	try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
	## Mehdi: add legend right on -- legend=levels(factor(factorLevels))
			
	try(dev.off())
    }
	return(retDf)	

}
