
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

kegg_pathway_GRS <- function(resultdf,species)
{

##	library(CLEAN)
	score <- 10^(-resultdf$geneTable[,2])
	org<-"Hs"
	if(species == "human") org<-"Hs"
	if(species == "mouse") org<-"Mm"
	if(species == "rat") org<-"Rn"
	res <- CLEAN::RandomSet(sigvals=score, geneids=resultdf$geneTable[,1], functionalCategories=c("KEGG"), species=org, minFDR=0.1, verbose=T, minGenesInCategory=10, maxGenesInCategory=5000)
	enrichPathways <- unique(as.character(res[[1]][,1]))

	sigIndex <- which(resultdf$geneTable[,2]> resultdf$EvalueNullDistrQ[2])
	allGenes <- paste(resultdf$geneTable[,"geneID"],collapse=",",sep="")

#	org<-"hsa"
#	if(species == "human") org<-"hsa"
#	if(species == "mouse") org<-"mmu"
#	if(species == "rat") org<-"rno"

##	library(RMySQL)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", host="10.165.4.231", port = 4040, password='public')
	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
	sql <- paste("SELECT Distinct LMD.* FROM GeneListMapping GLM,Category C,ListCategory LC,ListMetaData LMD
		Where	
		C.Name = \"KEGG\"
		AND LMD.Name in (",paste(paste("'",enrichPathways,"'",sep=""),collapse=",",sep=""),")
		AND LC.ListID = LMD.ID
		AND LC.CategoryID = C.ID
		AND GLM.ListID = LC.ListID
		AND GLM.GeneID in (",allGenes,")
		",sep="")
	
	rs <- DBI::dbSendQuery(mycon, sql)
	keggs <- DBI::fetch(rs,n=-1)
	if(nrow(keggs)==0) 
	{
		print("Pathways not found")
		return()
	}

	allGenes <- unlist(strsplit(allGenes,","))

	urlbegin <- "http://www.genome.jp/kegg-bin/mark_pathway_www?@"
	urlmid <- "/default%3dwhite/"
	urlend <- ""
	keggTable <- NULL

	significanceBarKegg <- rep("green",length(allGenes))
	significanceBarKegg[sigIndex] <- "red"
	table(significanceBarKegg)

	r <- NULL

	# For each pathway, find genes
	for(i in 1:dim(keggs)[1])
	{
		sql <- paste("select GeneID from GeneListMapping where ListID =",keggs[i,"ID"],sep="")	
		rs <- DBI::dbSendQuery(mycon, sql)
		g <- DBI::fetch(rs,n=-1)
	
		#kId<- sub("hsa","",keggs[i,"Name"])

		commonGenes <- intersect(g[,1], allGenes)
		commonIndex <- match(commonGenes, allGenes)
		significanceBartemp <- significanceBarKegg[commonIndex]
		sig <- length(which(significanceBartemp == "red"))
		
		url <- paste(urlbegin, keggs[i,"Name"],urlmid, sep = "")
		#url <- paste(urlbegin, org,kId,urlmid, sep = "")

		for (j in 1:length(commonGenes))
		{
			urltemp <- paste(commonGenes[j], "%09", significanceBartemp[j], "/", sep = "")
			urlend <- paste(urlend, urltemp, sep = "")
		}

		url <- paste(url, urlend, sep = "")
		urlend <- ""
		# Display urls.  or display the drop down list of Kegg pathways on results page
		displayUrl <- paste("<p><a href=",url,">", keggs[i,"Name"],"</a>",sep="")
                #count how many are significant (yellow)
		r<-rbind(r,c(displayUrl,keggs[i,"Description"],length(commonGenes),sig))



	}


	DBI::dbDisconnect(mycon)
	return(as.data.frame(r))
}
