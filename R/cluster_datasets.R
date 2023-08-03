
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

cluster_datasets <- function(sessionID, dataClust, nFactorLevels, ifCluster, currRowOrder, rowLable, rowsepindex,
			    colsepindex, rsd, pal, prop_distance, isPropNone, isCentered, path_to_write, 
			    prop, factorLevels, allColors, nColors, sigDB="ilincs_sigs", chost="ilincs.org", debuging=FALSE)
{
        colLabels <- NULL 
	isExtraColumn = FALSE 
	sepWidthVector <- as.vector(c(0.001,0.001))
	#create fake clustering objects for genes and samples. These will be used for treeView when clustering option is "none"
	rowClust <- fastcluster::hclust(dist(1:3))
	nGenes <- length(currRowOrder)
	rowClust$merge <- rbind(c(-1,-2),matrix(c(seq(-3,-nGenes),seq(1:(nGenes-2))),byrow=F,ncol=2))
	rowClust$height <- rep(1,nGenes-1)
	rowClust$order <- 1:nGenes

	colClust <- fastcluster::hclust(dist(1:3))
	nSamples <- dim(dataClust)[2] - 2
	colClust$merge <- rbind(c(-1,-2),matrix(c(seq(-3,-nSamples),seq(1:(nSamples-2))),byrow=F,ncol=2))
	colClust$height <- rep(1,nSamples-1)
	colClust$order <- 1:nSamples
## Mehdi save arguments for a test
#save(sessionID,dataClust,nFactorLevels, rowLable, rowsepindex, colsepindex, rsd, pal, prop_distance, ifCluster,
#isCentered,colClust,rowClust,currRowOrder,nGenes, sepWidthVector, path_to_write,isPropNone,
#prop, factorLevels, nFactorLevels, allColors, nColors,
#file=paste(path_to_write,"cluster_datasets.RData",sep=""))

	if(ifCluster == "both")
	{
		colClust<- fastcluster::hclust(dist(data.matrix(t(dataClust[currRowOrder,-(1:2)]))),method="complete")
		rowClust<- fastcluster::hclust(dist(data.matrix((dataClust[currRowOrder,-(1:2)]))),method="complete")
		colDendro <- as.dendrogram(colClust)
		rowDendro <- as.dendrogram(rowClust)
		create_treeview_files(sessionID,dataClust[currRowOrder,],ifCluster,isCentered,colClust,rowClust,path_to_write, sigDB=sigDB, chost=chost, debuging=debuging)
	} else if(ifCluster == "col") {
		colClust<- fastcluster::hclust(dist(data.matrix(t(dataClust[currRowOrder,-(1:2)]))),method="complete")
		colDendro <- as.dendrogram(colClust)
		rowDendro <- FALSE 
		create_treeview_files(sessionID,dataClust[currRowOrder,],ifCluster,isCentered,colClust,rowClust,path_to_write, sigDB=sigDB, chost=chost, debuging=debuging)
	} else if(ifCluster == "row"|| ifCluster == "none") {
		colDendro <- FALSE 
		rowDendro <- FALSE 
                if(!is.null(prop_distance))
                {
			workingDataClust <- dataClust[currRowOrder,-(1:2)]
			workingRowLables <- rowLable[currRowOrder]
			sampleNames <- unique(colnames(dataClust)[-(1:2)])
			isExtraColumn = TRUE 
			
			repeat_distance <- length(grep(prop_distance[1],prop_distance))
			if(repeat_distance == 1)  # there is only one sample. so do not stack genes.  Also show the grid lines
			{
				#show 0 as extra variable and add 50 to everything above zero 
			        pos <- which(prop_distance == 0)	
				prop_distance = as.character(prop_distance)
				colLabels <- c(prop_distance[1:pos-1],0,as.numeric(prop_distance[pos:length(prop_distance)])+50)
				rowsepindex <- as.vector(seq(1,dim(workingDataClust)[1],1))
				colsepindex <- as.vector(seq(1,length(colLabels),1))
			} else {  # not show the grid lines, only mark start and end distance and add column separation after each property.
				if(isPropNone) {
					NoOfSamples <- repeat_distance 
				} else {
					NoOfSamples <- nFactorLevels
				}
				partition <- length(prop_distance)/NoOfSamples
			        pos <- which(prop_distance == 0)	
				posPart <- which(prop_distance[1:partition] == 0)
				colLabelsPart= NULL
				tempcolLabelspart <- c(prop_distance[1:posPart[1]-1],0,as.numeric(prop_distance[posPart[1]:(partition/length(posPart))])+50)
				colLabelspart <- rep(tempcolLabelspart,length(posPart))

				colLabels <- rep(colLabelspart,NoOfSamples)
				minValue <- min(as.numeric(as.character(colLabels)))
				maxValue <- max(as.numeric(as.character(colLabels)))

				colLabArr <- c((min(as.numeric(as.character(colLabels)))+100),(max(as.numeric(as.character(colLabels)))-100),0)
				#go thru the array and keep only min and max and 0. make everything else ' '
				for (c in 1:length(colLabels))
				{
					if(!is.element(as.numeric(as.character(colLabels[c])),colLabArr))
					{
						colLabels[c] <- ' '
					}
				}
				colLabels[which(colLabels == (minValue+100))] <- minValue
				colLabels[which(colLabels == (maxValue-100))] <- maxValue

				partition <- partition + length(posPart) 
				colsepindex <- seq(partition,length(colLabels), by= partition) 
				rowsepindex <- as.vector(seq(1,dim(workingDataClust)[1],1))
				
				sepWidthVector <- c(1,0.05)
			}
                } else { # no distance vector
                        colLabels <- NA
			sepWidthVector <- c(0.005,0.05)
                }
		if(isExtraColumn)
		{	
			t1 <- Sys.time()
			insert <- rep(NA,nGenes)
			newClust <- as.vector(as.matrix(dataClust[currRowOrder,1:2])) 
			tempClust <- as.vector(as.matrix(dataClust[currRowOrder,-(1:2)]))	
			colnam <- colnames(dataClust)[-(1:2)]
			start <- 1
			posOriginal <- pos
			pos <- (pos-1)*nGenes
			for( i in 1:(length(pos)))
			{
				newClust <- c(newClust,tempClust[start:pos[i]],insert)
				start <- pos[i]+1
				rsd <- c(rsd[1:posOriginal[i]-1],rsd[posOriginal[i]-1],rsd[posOriginal[i]:length(rsd)])
				colnam <- c(colnam[1:posOriginal[i]-1],"Separator",colnam[posOriginal[i]:length(colnam)])
			}
			newClust <- c(newClust,tempClust[(pos[i]+1):length(tempClust)])
			dataClust <- matrix(newClust,nGenes,(ncol(dataClust)+length(pos)))
			#dataClust <- cbind(dataClust[,1:2],as.numeric(dataClust[,3:ncol(dataClust)]))
			dataClust <- as.data.frame(dataClust)
			#for(i in 3:ncol(dataClust)) dataClust[,i] <- as.numeric(as.character(dataClust[,i]))
			z <- dataClust[,3:ncol(dataClust)]
			z <- apply(z,2,function(x) x <- as.numeric(x))
			dataClust <- cbind(dataClust[,1:2],z)

			colnam <- c(colnames(dataClust)[1:2],colnam)
			colnames(dataClust) <- colnam

			#create fake clustering object for samples.
			colClust <- fastcluster::hclust(dist(1:3))
			nSamples <- dim(dataClust)[2] - 2
			colClust$merge <- rbind(c(-1,-2),matrix(c(seq(-3,-nSamples),seq(1:(nSamples-2))),byrow=F,ncol=2))
			colClust$height <- rep(1,nSamples-1)
			colClust$order <- 1:nSamples
			t2 <- Sys.time()
		
		}#end ifExtraColumn
		if(ifCluster == "row")
		{
			rowClust <- fastcluster::hclust(dist(data.matrix((dataClust[currRowOrder,-(1:2)]))),method="complete")
			rowDendro <- as.dendrogram(rowClust)
		}
		
		create_treeview_files(sessionID,dataClust[currRowOrder,],ifCluster,isCentered,colClust,rowClust,path_to_write, sigDB=sigDB, chost=chost, debuging=debuging)
	} # end of if row or none 

	if(!is.null(prop_distance))
	{
	    if(isCentered) title <- "Centered Promoter Profile" else title <- "Original Promoter Profile"
	} else {
	     if(isCentered) title <- "Clustering centered log2 levels" else title <- "Clustering absolute log2 levels" 
	}
	rowcolFont <- min(0.2 + 1/log10(length(currRowOrder)), 0.7)#1)
#	RLmargins <- max(nchar(rowLable))
	RLmargins <- min(50, (max(nchar(rowLable)) + max(nchar(levels(factor(factorLevels))))) * 0.7)
	if (isPropNone) RLmargins = max(nchar(rowLable)) * 0.7
	TBmargins <- min(25, max(nchar(colnames(dataClust))) * 0.7) # max(nchar(colLabels))
	#save(dataClust, rowLable, colLabels, currRowOrder, file=paste(path_to_write,"/datacluster.RData",sep=""))
	#print(paste("*********", TBmargins, RLmargins))

    servSet <- getServerSettings(debuging=debuging, chost=chost)
    if(servSet$compute == "dev.ilincs.org") {
	if(isPropNone)
	{
		if(isCentered)
		{
			if(!is.null(prop_distance))  # colsepindex only for ChIP datatype
			{
			try(gimmHeat(data.matrix(dataClust[currRowOrder,-(1:2)]),gimmTrunc=c(-3,3),trace="none",tracecol = "none",density="none",colsep = colsepindex,sepwidth=sepWidthVector,sepcol="grey40",col = pal,labCol=colLabels,labRow=rowLable[currRowOrder],scale="none",main=title, margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))# changed 0.6 to 1.0 in cexRow
#			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				##pch=rep(19,nFactorLevels), ncol=3, bty="o", bg=NA, pt.cex=2)## Mehdi added legend margins=c(12,12) (8,15)
#			legend(x=0.1,y=1,title = prop,legend=factorLevels,col=allColors[(((1:nFactorLevels)-1)%%nColors)+1],
#				pch=rep(19,nFactorLevels),ncol=3,bty="o",bg = "lightgrey",pt.cex=2)
			} else {
			try(gimmHeat(data.matrix(dataClust[currRowOrder,-(1:2)]),gimmTrunc=c(-3,3),trace="none",tracecol = "none",density="none",col = pal,labCol=colLabels,labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
#			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				## Mehdi
			}
		} else {
			if(!is.null(prop_distance))  # colsepindex only for ChIP datatype
			{#print(paste("************************this is it************************"))

			try(gplots::heatmap.2(data.matrix(dataClust[currRowOrder,-(1:2)]),trace="none",density="none",colsep = colsepindex,sepwidth=sepWidthVector,sepcol="grey40",col = pal,labCol=colLabels,labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
# 			print(paste(colLabels))
# 			print(paste(rowLable))
# 			colLabels<-colLabels
				## Mehdi
			} else {
			try(gplots::heatmap.2(data.matrix(dataClust[currRowOrder,-(1:2)]),trace="none",density="none",col = pal,labCol=colLabels,labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
#			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				## Mehdi
			}
		}
	} else {
		if(isCentered)
		{
			if(!is.null(prop_distance))
			{
			try(gimmHeat(data.matrix(dataClust[currRowOrder,-(1:2)]),gimmTrunc=c(-3,3),trace="none",tracecol = "none",density="none",colsep = colsepindex,sepwidth=sepWidthVector,sepcol="grey40",col = pal,ColSideColors=rsd,labCol=colLabels, labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				##pch=rep(19,nFactorLevels), ncol=3, bty="o", bg=NA, pt.cex=2)## Mehdi legend
#			legend(x=0.1,y=1,title = prop,legend=factorLevels,col=allColors[(((1:nFactorLevels)-1)%%nColors)+1],
#				pch=rep(19,nFactorLevels),ncol=3,bty="o",bg = "lightgrey",pt.cex=2)
			} else {
			try(gimmHeat(data.matrix(dataClust[currRowOrder,-(1:2)]),gimmTrunc=c(-3,3),trace="none",tracecol = "none",density="none",col = pal,ColSideColors=rsd,labCol=colLabels, labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
			}
			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				## Mehdi
		} else {
			if(!is.null(prop_distance))
			{
			try(gplots::heatmap.2(data.matrix(dataClust[currRowOrder,-(1:2)]),trace="none",density="none",colsep = colsepindex,sepwidth=sepWidthVector,sepcol="grey40",col = pal,ColSideColors=rsd,labCol=colLabels,labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				## Mehdi
			} else {
			try(gplots::heatmap.2(data.matrix(dataClust[currRowOrder,-(1:2)]),trace="none",density="none",col = pal,ColSideColors=rsd,labCol=colLabels,labRow=rowLable[currRowOrder],scale="none",main=title,margins= c(TBmargins,RLmargins), dendrogram=ifCluster,Colv=colDendro,Rowv=rowDendro, keysize = 0.8, lhei = c(1,max(2,(0.05 * length(currRowOrder)))), cexCol = rowcolFont, cexRow = rowcolFont))
			try(legend("topright",title = prop, legend=factorLevels, fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], cex=0.7))
				## Mehdi  fill=allColors[(((1:nFactorLevels)-1)%%nColors)+1],col=allColors[(((1:nFactorLevels)-1)%%nColors)+1], 
				#legend=levels(factor(factorLevels))
				#print(paste("###############", factorLevels, "##############"))
			}
		}
	}
    }
}
