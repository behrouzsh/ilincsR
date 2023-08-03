
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

plot_nbr <- function(db,exp,cdtdata,prop,web=F, debuging=FALSE){
	pal <- marray::maPalette(low="blue",high="yellow",mid="black")

	meas_limits <- initiate(db,exp)
	sampledata <- get_sample_properties(db=db,exp=exp,meas_limits=meas_limits, debuging)
	dataTable <- parse_sample_properties(propertyTags=sampledata)
        propertydata <- data.frame(dataTable[,1],dataTable[,prop])
        graphDir <- "../htdocs/tmp/"
        graphURLroot <- "/tmp/"
        allColors <- c("blue","red","green","yellow","purple","brown","black","blue4","red4","green4","yellow4")
        matchingSamples <- match(as.character(unlist(propertydata[,1])),unlist(colnames(cdtData[,-(1:4)])))

        factorLevels <- sort(unique(propertydata[,2]))
	nFactorLevels <- length(factorLevels)

        nObs <- dim(cdtData)[1]
        nCols <- dim(cdtData)[2]
        rsd <- rep("blue",nCols-4)
                for(i in 1:nFactorLevels){
                        rsd[propertydata[,2]==factorLevels[i]] <- allColors[i]
                }

                ss <- apply(cdtData,2,function(x) as.numeric(x))
                colDist <- as.dist((1-cor(data.matrix(ss[,5:nCols])))/2)
                colClust <- hclust(colDist,method="average")
                rowDist <- as.dist((1-cor(t(data.matrix(ss[,5:nCols]))))/2)
                rowClust <- hclust(rowDist,method="average")

         if(web) webPNG("temp.png",res=120,height=6,width=6,typ="png16m")
         gplots::heatmap.2(ss[,5:nCols],trace="none",Colv=as.dendrogram(colClust),Rowv=as.dendrogram(rowClust),ColSideColors=rsd,col=pal,labCol=NA,labRow=as.character(cdtData[,3]),margins=c(8,15)) ## (5,10)
         if(web) img(src = "temp.png")
	plotSize <- floor(nFactorLevels/2)
         if(web) webPNG("temp2.png",res=120,height=plotSize,width=6,typ="png16m")
        par(mar=c(0,5,3,5))
          plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",main=paste("Variable ",prop,sep=""),bty="n",cex.main=2)
          legend(x=0.3,y=1,legend=factorLevels,col=allColors[1:nFactorLevels],pch=rep(19,nFactorLevels),ncol=1,bty="n",cex=2)
         if(web) img(src = "temp2.png")
}
