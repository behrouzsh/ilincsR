
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

create_expr_for_gff_sample <- function(result,up,down,window.size=100,missingWindow=0)
{	
##	require(IRanges)
##	require(Biobase)

	# recreate the genome
	genome = unique(result[,1:8])
	allChrs = as.character(na.omit(unique(genome$chrom)))
	n.scores = ceiling((up+down)/window.size)
	n.refseqs = length(unique(genome$RefID))
	exps = matrix(missingWindow,nrow=n.refseqs,ncol=n.scores)
	rownames(exps) = unique(genome$RefID)
	exps = t(exps)
	result = result[order(result$chrom,result$chromStart),]
	

	for(chr in allChrs)
	{
		xset <- result[result$chrom == chr,]
		genome.0 <- which(genome$chrom == chr)
		window.starts <- as.vector(sapply(genome.0,function(x)
		{
			if(genome$Strand[x]=="+") {
				seq(from=genome$TSS[x]-up,to=genome$TSS[x]+down-1,by=window.size)
			 } else {
				seq(from=genome$TSS[x]+up-window.size,to=genome$TSS[x]-down,by=-window.size)}
		}))
		window.ends <- window.starts + window.size
	
		window.pos <- IRanges::IRanges(window.starts,window.ends)
		probe.pos <- IRanges::IRanges(xset$chromStart,xset$chromEnd)
		window.scores <- rep(missingWindow,length(window.starts))
		if(length(window.pos)>0 & length(probe.pos)>0)
		{
#@			tree <- IRanges::IntervalTree(window.pos)
#@@			tree <- NCList(IRanges(genome.0.starts,genome.0.ends))
			tree <- IRanges::NCList(window.pos)
#@			tree.overlap <- as.matrix(overlap(tree,probe.pos))
			tree.overlap <- as.matrix(IRanges::findOverlaps(tree,probe.pos))
##			tables <- as.matrix(findOverlaps(query=IRanges(temp$Start,temp$End),tree,type="within",select="all", algorithm="nclist"))
			if(dim(tree.overlap)[1]>0)
			{
				matchings <- as.data.frame(tree.overlap)
				scores <- tapply(matchings$query,matchings$subject,function(x)mean(xset$dataValue[x],na.rm=T))
				window.scores[as.numeric(names(scores))]=scores
			}
		}
		exps[1:n.scores,genome.0] = as.vector(window.scores)
	}
	exps[is.nan(exps)] <- missingWindow
	exps <- t(exps)
	exps <- cbind(exps,as.vector(rep(result[1,"Sample"],dim(exps)[1])))
	#result.eset = new("ExpressionSet",exprs=exps)
	#featureNames(result.eset)=genome$RefID
	#fData(result.eset)=data.frame(genome[,1:2],row.names=featureNames(result.eset))
	#sampleNames(result.eset)=seq(-up,down,by=window.size)[1:n.scores]
	#pData(result.eset)=data.frame(distance=sampleNames(result.eset),row.names=sampleNames(result.eset))
	#result.eset
	return(exps)
}
