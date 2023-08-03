
#' A function for ...
#'
#' This function allows you to find Human homolog genes for non-Human specious.
#' @param glist A vector of geneIDs.
#' @param tax TaxID for the specious ...
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

# # getHomologousGenes <- function(glist){
# # 	tax <- "9606"
# # #	glist<-unlist(strsplit(glist,","))
# # #	glist<-sub(" ","",glist)
# # #	glist<-gsub(",,",",",glist)
# #  #       glist<-gsub("'NA',","",glist)
# #  #       glist<-gsub(",'NA'","",glist)
# # 
# # # 	load(url("http://www.eh3.uc.edu/geneinfo.RData"))
# # 	data(geneinfo)
# # 	if(sum(as.numeric(glist), na.rm=T) == 0) {
# # 	    glist <- geneinfo[match(glist, geneinfo$Symbol),"GeneID"]
# # 	}
# # 	genedf <- geneinfo[match(glist,geneinfo[,1]),]
# # 	l <- which(genedf[,"TaxID"] == tax)
# # 	## homologenes
# # 	hml <- genedf[genedf["TaxID"]!=tax,"HomologeneID"]
# # 	x <- geneinfo[geneinfo["TaxID"]==tax,]
# # 	ind <- match(hml,x[,"HomologeneID"])
# # 	hgenedf <- x[ind,]
# # 	if(length(l) > 0) {
# # 		hgenedf[l,] <- genedf[l,]
# # 	}
# # 	retArr <- NULL
# # 	retArr <- cbind(genedf[,"GeneID"], hgenedf[,"GeneID"])
# # 	colnames(retArr) <- c("geneID","homoloGeneID")
# # 	return(retArr)
# # }


# # # # # getHomologousGenes <- function(glist, org = "Hs"){
# # # # # 	
# # # # # 
# # # # # 	tax <- org2tax(org)
# # # # # 	data(geneinfo)
# # # # # 	if(length(na.omit(as.integer(glist))) < (length(glist)/2)) {
# # # # # 	    glist <- geneinfo[match(glist, geneinfo$Symbol),"GeneID"]
# # # # # 	}
# # # # # 	glist <- glist[!is.na(glist)]
# # # # # 	glist <- glist[glist %in% geneinfo$GeneID]
# # # # # 	
# # # # # # 	originals <- geneinfo[geneinfo$GeneID %in% glist, ]
# # # # # # 	homoIDs <- unique(originals$HomologeneID)
# # # # # # 	res <- geneinfo[geneinfo$HomologeneID %in% homoIDs & as.character(geneinfo$TaxID) == tax,]
# # # # # # 	retArr <- res[,c("HomologeneID", "GeneID")]
# # # # # # 	colnames(retArr) <- c("GeneID", "HomologeneID")
# # # # # 	
# # # # # 	genedf <- geneinfo[match(glist,geneinfo$GeneID),]
# # # # # 	l <- which(genedf["TaxID"] == tax)
# # # # # 	## homologenes
# # # # # 	hml <- genedf[genedf["TaxID"]!=tax,"HomologeneID"]
# # # # # 	x <- geneinfo[geneinfo["TaxID"]==tax,]
# # # # # 	ind <- match(hml,x[,"HomologeneID"])
# # # # # 	hgenedf <- x[ind,]
# # # # # 	if(length(l) > 0) {
# # # # # 		hgenedf[l,] <- genedf[l,]
# # # # # 	}
# # # # # 	retArr <- NULL
# # # # # 	retArr <- cbind.data.frame(genedf[,"GeneID"], hgenedf[,c("GeneID", "Symbol")])
# # # # # 	colnames(retArr) <- c("geneID","homoloGeneID","homoloGeneSymbol")
# # # # # 	return(retArr)
# # # # # }


# servSet <- getServerSettings(debuging=T)
# mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GeneDB", port = servSet$port, host = servSet$host, password = "public")
# geneinfo2 <- DBI::dbGetQuery(mycon, "SELECT * FROM GeneInfo")


getHomologousGenes <- function(glist, org = "Hs"){
	

	tax <- org2tax(org)
	data(geneinfo)
	if(length(na.omit(as.integer(glist))) < (length(glist)/2)) {
	    glist <- geneinfo[match(glist, geneinfo$Symbol),"GeneID"]
	}
	glist <- glist[!is.na(glist)]
	glist <- glist[glist %in% geneinfo$GeneID]
	
# 	originals <- geneinfo[geneinfo$GeneID %in% glist, ]
# 	homoIDs <- unique(originals$HomologeneID)
# 	res <- geneinfo[geneinfo$HomologeneID %in% homoIDs & as.character(geneinfo$TaxID) == tax,]
# 	retArr <- res[,c("HomologeneID", "GeneID")]
# 	colnames(retArr) <- c("GeneID", "HomologeneID")
	
	genedf <- geneinfo[match(glist,geneinfo$GeneID),]
# 	l <- which(genedf["TaxID"] == tax)
	l <- which(genedf["TaxID"] != tax)
	## homologenes
# 	for(j in l) {
# 	    genedf[l,] <- 
# 	}
	
	hml <- genedf[genedf["TaxID"]!=tax,"HomologeneID"] # na.omit()
	x <- geneinfo[geneinfo["TaxID"]==tax,]
	ind <- match(hml,x[,"HomologeneID"])
	hgenedf <- genedf
	hgenedf[l,] <- x[ind,]
# 	if(length(l) > 0) {
# 		hgenedf[l,] <- genedf[l,]
# 	}
# # # 	if(length(hml) > 0) {
# # # 		hgenedf[l,] <- genedf[l,]
# # # 	}
	retArr <- NULL
# 	retArr <- cbind.data.frame(genedf[,"GeneID"], hgenedf[,c("GeneID", "Symbol")])
	retArr <- cbind.data.frame(genedf[,c("GeneID", "Symbol")], hgenedf[,c("GeneID", "Symbol")])
	colnames(retArr) <- c("geneID","Symbol","homoloGeneID","homoloGeneSymbol")
	return(retArr)
}

