
#' A function for ...
#'
#' This function allows you to ...
#' @param sigFile A signature file from ilincs.
#' @param libName The library user wants to find enrichments from (currently LIB_5 for drugCategories and LIB_6 for genesCategories).
#' @param userUploadedProfilePath Path to the signature file.
#' @param sigDB The database for signature libraries ("sigs or new").
#' @param two.sided How to calculate p-values for Zscores.
#' @param authors M. Fazel-Najafabadi
#' @keywords mehdi_wtcor
#' @export
#' @examples
#' res <- sigEnrichment(sigFile="sig_Fri_Sep_28_16_38_27_2018_2976102.xls", libName="LIB_6", debuging=T)
#'


sigEnrichment <- function(sigFile, libName, userUploadedProfilePath="/mnt/raid/tmp/", logPcut=10, sigDB="ilincs_sigs", two.sided=FALSE, org="Hs", debuging=FALSE) {
    res <- list(Remark="Done", enrichment=NA)
    if (!(libName %in% c("LIB_5", "LIB_6"))) {
	res$Remark <- "Selected signature library is not appropriate!"
	return(res)
    }
    if (length(grep("processedSig_", sigFile)) == 1) {
      queryProfile <- try(read.table(file = paste0(userUploadedProfilePath,"/", sigFile), header=TRUE, sep="\t", stringsAsFactors=FALSE)) # , quote=""
    } else {
      queryProfile <- uploadedSigProcess(sigFile, userUploadedProfilePath, write=FALSE, debuging=debuging)$Sig
    }
    if(dim(queryProfile)[2] == 3) queryProfile$Significance_pvalue <- 1
    colnames(queryProfile) <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
    queryProfile <- queryProfile[, c("ID_geneid", "Value_LogDiffExp", "Significance_pvalue")]
    homoloArr <- getHomologousGenes(as.character(queryProfile[,1]), org=org)
    if(length(!is.na(homoloArr[,3])) == 0){
	    res$Remark <- "No homolo genes"
    }
    queryProfile[,1] <- homoloArr[match(queryProfile[,1], homoloArr[,1]),3]
    queryProfile <- queryProfile[which(!is.na(queryProfile[,1])),]
    if(length(unique(queryProfile[,"ID_geneid"])) < dim(queryProfile)[1]) {
	if(dim(queryProfile)[2] == 3) unqIdPvals <- split(queryProfile[,"Significance_pvalue"], queryProfile[,"ID_geneid"])
	if(dim(queryProfile)[2] == 3) mergePvals <- unlist(lapply(unqIdPvals, function(x) exp(mean(log(x), na.rm = TRUE)))) ## Mehdi: I changed it so we dont loose genes with p-val=NA
	unqIdCoeffs <- split(queryProfile[,"Value_LogDiffExp"], queryProfile[,"ID_geneid"])
	mergeCoeffs <- unlist(lapply(unqIdCoeffs, function(x) mean(as.numeric(x), na.rm = TRUE)))
	if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(ID_geneid=names(unqIdCoeffs), Value_LogDiffExp=mergeCoeffs, Significance_pvalue=mergePvals, stringsAsFactors=FALSE)
	if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(ID_geneid=names(unqIdCoeffs), Value_LogDiffExp=mergeCoeffs, Significance_pvalue=1, stringsAsFactors=FALSE)
    }
    rownames(queryProfile) <- queryProfile$ID_geneid
    a <- Sys.time()
    servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
    mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=servSet$sigDB, host=servSet$host, port = servSet$port, password="public")
    sql <- paste0("select diffexpTablePath,pvaluesTablePath from signatureLibraries where libraryID = '",libName,"'")
    rs <- DBI::dbGetQuery(mycon, sql)
    libraryPath <- rs[1,"diffexpTablePath"]
    libraryProb <- rs[1,"pvaluesTablePath"]
    b <- Sys.time()
    refProfiles <- get(load(libraryPath))
    genes <- intersect(rownames(queryProfile), rownames(refProfiles))
    queryData <- queryProfile[genes, ]
    refProfiles <- refProfiles[genes, ]
    if(dim(refProfiles)[1] < 3) {
	res$Remark <- "Very few genes!"
	return(res)
    }
    libraryProb <- sub("2018/", "2018/signed_pLogs/", libraryProb)
    libraryProb <- sub("PValues", "Plogs", libraryProb)
    refProbs <- get(load(libraryProb))
    refProbs <- refProbs[genes, ]
    queryData[is.na(queryData$Significance_pvalue),"Significance_pvalue"] <- 1
    log10p <- -log10(queryData$Significance_pvalue)
    log10p[log10p > logPcut] <- logPcut
    wt <- abs(refProbs) + log10p
    AllCor <- mehdi_wtcor(queryData$Value_LogDiffExp, refProfiles, wt, topSigs=FALSE)
    AllCor <- data.frame(signatureID=colnames(AllCor), t(AllCor)[,c(1,4,3)], stringsAsFactors=F)
    AllCor[,4] <- AllCor[,4] + 2
    AllCor[,3][AllCor[,3]==0] <- 2.225074e-308
    colnames(AllCor) <- c("signatureID", "Similarity", "PValue", "nGene")
    AllCor$p_up <- ifelse(AllCor[, 2] >= 0, AllCor[, 3]/2, 1 - (AllCor[, 3]/2))
    AllCor$p_dn <- ifelse(AllCor[, 2] <= 0, AllCor[, 3]/2, 1 - (AllCor[, 3]/2))
    AllCor$p_up[AllCor$p_up == 0] <- min(AllCor$p_up[AllCor$p_up != 0]) * 0.01
    AllCor$p_dn[AllCor$p_dn == 0] <- min(AllCor$p_dn[AllCor$p_dn != 0]) * 0.01
    tmpCor <- if(libName=="LIB_5") drugCategory$meta else geneCategory$meta
    mch <- match(AllCor$signatureID, tmpCor$signatureID)
    AllCor$category <- tmpCor[mch, "category"]
    AllCor$compound <- tmpCor[mch, "compound"]
    AllCor$p_up <- -log10(AllCor$p_up)
    AllCor$p_dn <- -log10(AllCor$p_dn)
    AllCor0 <- AllCor[!is.na(AllCor$p_up),]
    AllCorup <- split(AllCor0$p_up, AllCor0$category)
    sp <- unique(tmpCor[,5:6])
    rownames(sp) <- sp$category
    tmpup <- data.frame(categoryName=names(AllCorup), compound=sp[names(AllCorup), "compound"], 
			m=sapply(AllCorup, length), G=dim(AllCor0)[1], xBar=sapply(AllCorup, mean), mu=mean(AllCor0$p_up), stringsAsFactors=FALSE)
    scores2 <- sum(AllCor0$p_up^2, na.rm = TRUE)
    tmpup$sigma <- sqrt((scores2/tmpup$G - tmpup$mu^2) * (tmpup$G - tmpup$m)/(tmpup$G - 1)/tmpup$m)
    tmpup$zScore <- (tmpup$xBar - tmpup$mu)/tmpup$sigma
    tmpup$connDirection <- "+"
    AllCor0 <- AllCor[!is.na(AllCor$p_dn),]
    AllCordn <- split(AllCor0$p_dn, AllCor0$category)
    sp <- unique(tmpCor[,5:6])
    rownames(sp) <- sp$category
    tmpdn <- data.frame(categoryName=names(AllCordn), compound=sp[names(AllCorup), "compound"], 
			m=sapply(AllCordn, length), G=dim(AllCor0)[1], xBar=sapply(AllCordn, mean), mu=mean(AllCor0$p_up), stringsAsFactors=FALSE)
    scores2 <- sum(AllCor0$p_up^2, na.rm = TRUE)
    tmpdn$sigma <- sqrt((scores2/tmpdn$G - tmpdn$mu^2) * (tmpdn$G - tmpdn$m)/(tmpdn$G - 1)/tmpdn$m)
    tmpdn$zScore <- (tmpdn$xBar - tmpdn$mu)/tmpdn$sigma
    tmpdn$connDirection <- "-"
    enrichment <- rbind(tmpup, tmpdn)
    if(two.sided) {
	    enrichment$pValue <- pnorm(-abs(enrichment$zScore)) * 2
    } else {
	    enrichment$pValue <- pnorm(-enrichment$zScore)
    }
    enrichment <- enrichment[,c("categoryName","compound", "zScore", "pValue", "m", "connDirection")]
    colnames(enrichment)[5] <- "nSignature"
    enrichment$FDR <- p.adjust(enrichment$pValue, "fdr")
    enrichment <- enrichment[,c(1,2,3,4,7,5,6)]
    enrichment <- enrichment[enrichment$FDR <= 0.1,]
    
    if(dim(enrichment)[1] < 1) {
	res$Remark <- "No enriched signatures found!"
	##return(res)
    }
    enrichment <- enrichment[order(abs(enrichment$zScore),decreasing=T),]
    rownames(enrichment) <- NULL
    if(libName == "LIB_5") {
	perturbagenIDs <- paste0(enrichment$categoryName, collapse="','")
	rs <- DBI::dbGetQuery(mycon, paste0("select distinct perturbagenID,lincsPertID from signaturesMeta3 where perturbagenID in ('", perturbagenIDs, "')"))
	enrichment$lincsPertID <- rs[na.omit(match(enrichment$categoryName, rs$perturbagenID)),"lincsPertID"]
	targets <- drugTargets(drugID=enrichment$categoryName, sigDB=sigDB, debuging=debuging)
	enrichment$Targets <- targets$targets[match(enrichment$categoryName, targets$drugID)]
	enrichment <- enrichment[, c("categoryName", "compound", "lincsPertID", "Targets", "zScore", "pValue", "FDR", "nSignature", "connDirection")]
    }
    if(libName == "LIB_6") {
        if(length(na.omit(as.integer(enrichment[, 1]))) < length(enrichment[, 1])/2) {
	    data(geneinfo)
	    ids <- as.character(geneinfo[match(enrichment[, 1], geneinfo$Symbol),"GeneID"])
	}
	keggPaths <- get(data(list=paste0("keggPathways", org)))
	pathwayidlist <- gsub(" - Homo sapiens [(]human[)]", "",get(data("keggPathways"))[[org]])
	enrichment$pathways <- sapply(ids, function(j) {
	    if(!is.na(j)) {
		genesIn <- lapply(keggPaths, function(k) k[[1]]$GENE[k[[1]]$GENE %in% j])
		paste(pathwayidlist[names(unlist(genesIn))], collapse=" | ")
	    } else NA
	})
    }
    DBI::dbDisconnect(mycon)
    res$enrichment <- enrichment    
    res
}

