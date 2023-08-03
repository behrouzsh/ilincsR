
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

kegg_pathways <- function(gpgene,gpprobe,dataclust2,pkeggvalue)
{
##		library(RMySQL)

                keggGeneID <- as.vector(unlist(strsplit(gpgene,"DELIM")))
		s<-NULL
		for(a in keggGeneID)
		{
			if(is.null(s))
			{
				s<-paste("'",a,"'",sep="")
			} else {
				s<-paste(s,",'",a,"'",sep="")
			}

		}
		#sql<-paste("SELECT distinct GD.Name,GD.Description,GD.List FROM GeneDB.GeneToListMap GM,GeneDistance.KeggLists GD WHERE GM.GeneID in (",s,") and GM.GeneListTable = 'KeggLists' and GM.GeneList = GD.Name",sep="")
		sql<-"Select GLM.GeneID,LMD.Name,LMD.Description from ListMetaData LMD, GeneListMapping GLM, ListCategory LC, Category C where C.Name = \"KEGG\" And LC.CategoryID = C.ID And LMD.ID = LC.ListID AND GLM.ListID = LC.ListID"
                mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="GeneDB", port = 4040,host="10.165.4.231", password='public')
                rs <- DBI::dbSendQuery(mycon, sql)
                g <- DBI::fetch(rs,n=-1)
                DBI::dbDisconnect(mycon)
		if(nrow(g)==0) 
		{

			return()
		} else {
		
		q <- split(g[,1],g[2])
		y <- lapply(q, function(x) paste(unlist(x),collapse=":",sep=""))
		dd <- split(g[,3],g[,2])

		desc <- lapply(dd,function(x) x[1])
		keggs <- data.frame(I(names(desc)),I(desc),I(y) )
		keggs <- apply(keggs,c(1,2),function(x) x <- as.character(x))
		colnames(keggs)<-c("Name","Description","List")
	
                pathwayArray <- keggs[,1]   #kegg pathway name
                pathwayGenelist <- keggs[,3] #corresponding gene list

#               Construct url for all pathways according to color generated from Pvalue (for each gene)
                dataClustprobe <- apply(as.matrix(unlist(strsplit(dataclust2,"DELIM"))), 1, function(x) paste(unlist(strsplit(x,":"))[1], sep = ":"))
                keggProbe <- as.vector(unlist(strsplit(gpprobe,"DELIM")))
##                require(gdata)
                keggIndex <- match(gdata::trim(dataClustprobe), gdata::trim(keggProbe))
		Pvalues<-as.double(unlist(strsplit(as.character(pkeggvalue),"DELIM")))
                PvaluesKegg <- (Pvalues[keggIndex])

                significanceBarKegg <- rep("green", length(Pvalues))
                significanceBarKegg[PvaluesKegg<0.05]<-"red"
                table(significanceBarKegg)

                urlbegin <- "http://www.genome.jp/kegg-bin/mark_pathway_www?@"
                urlmid <- "/reference%3dwhite/"
                urlend <- ""
                keggTable <- NULL
                for(i in 1:length(pathwayArray))
                {
                        url <- paste(urlbegin, pathwayArray[i],urlmid, sep = "")
                        # get intersection of genes in kegg pathway and query gene list. Contruct url only using those genes
                        # gene list in kegg pathway is ":" separated
                        keggpathwaylist <- unlist(strsplit(pathwayGenelist[i],":"))
                        commonGenes <- intersect(keggpathwaylist, keggGeneID)
                        commonIndex <- match(commonGenes, keggGeneID)
                        significanceBartemp <-significanceBarKegg[commonIndex]

                        for (j in 1:length(commonGenes))
                        {
                                urltemp <- paste(commonGenes[j], "%09", significanceBartemp[j], "/", sep = "")
                                urlend <- paste(urlend, urltemp, sep = "")
                        }
                        url <- paste(url, urlend, sep = "")
                        urlend <- ""
#                       Display urls.  or display the drop down list of Kegg pathways on results page
                        displayUrl <- paste("<p><a href=",url,">", pathwayArray[i],"</a>",sep="")
                        #count how many are significant (yellow)
                        keggSigGenes <- length(which(significanceBartemp == "red"))
                        keggTable <- rbind(keggTable, c(displayUrl,keggs[i,2],length(keggpathwaylist), length(commonGenes), keggSigGenes))
                }
                colnames(keggTable) <- c("Kegg pathway","Description", "Total number of genes in the Kegg pathway"," # query genes in the pathway", " # significant Genes in the pathway")

		kt <- keggTable[order(as.numeric(keggTable[,5]), decreasing=TRUE),]
		return(as.data.frame(kt))
		}

}
