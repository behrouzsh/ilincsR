
#' A function to preProcess an uploaded signature file
#'
#' This function allows you to translate gene IDs to gene Names and vice versa and remove unavailable genes in iLincs portal.
#' @param exp Your uploaded file name. The file format should be in plain text with tab/comma/semicolon/space delimiters (only one of them in a file though)
#' @param userUploadedProfilePath The folder path containing your file.
#' @keywords Entrez Gene ID
#' @export 
#' @examples
#' processedSig <- uploadedSigProcess(exp = "sample.xls", userUploadedProfilePath="/mnt/raid/tmp/")


uploadedSigProcess <- function(exp=NULL, userUploadedProfilePath="/mnt/raid/tmp/", path_to_write="/mnt/raid/tmp/", org="Hs", write=TRUE, topGenes=NULL, debuging=FALSE) {
		
#@!    if(!is.null(userUploadedProfilePath)){
	    #upload custom profile from a file
	    res <- list(sessionID = NA, 
			orgFileName = exp, 
			fileName = NA, 
			shortfileName = NA, 
			Remark = "Error, wrong format or missing field", 
			Pval = TRUE, 
			entries = NA, 
			notFound1 = NA,
			notFound2 = NA,
			foundGenes = NA, 
			foundUniqueGenes = NA, 
			Sig = NA, 
# 			shortSig = NA, 
			commonGenes = NA#,
# 			shortcommonGenes = NA
			)
			
	    if(exp=="UpDn") {
		upgenes <- topGenes$genesUp # [[1]]
		dngenes <- topGenes$genesDown # [[2]]
		res$notFound2 <- c(upgenes, dngenes)
		
		wrongs <- intersect(upgenes, dngenes)
		if(length(wrongs)>0) {
		    res$Remark <- paste0(res$Remark, " - Gene(s) ", paste(wrongs, collapse="::"), " caused confusion!")
		    upgenes <- upgenes[!(upgenes %in% wrongs)]
		    dngenes <- dngenes[!(dngenes %in% wrongs)]
		}
		data(geneinfo)
		if(sum(as.numeric(upgenes), na.rm=T) == 0) upgenes <- geneinfo[match(upgenes, geneinfo$Symbol),"GeneID"]
		if(sum(as.numeric(dngenes), na.rm=T) == 0) dngenes <- geneinfo[match(dngenes, geneinfo$Symbol),"GeneID"]
		queryProfileUP <- data.frame(geneID=upgenes, coefficients=1, stringsAsFactors=FALSE)
		queryProfileDN <- data.frame(geneID=dngenes, coefficients=-1, stringsAsFactors=FALSE)
		queryProfile <- rbind.data.frame(queryProfileUP, queryProfileDN)
	    } else {
		queryProfile <- list()
		queryProfile[[1]] <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE , sep="\t", stringsAsFactors=FALSE), silent=TRUE) # , quote=""
		queryProfile[[2]] <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, sep=",", stringsAsFactors=FALSE), silent=TRUE) # , quote=""
		queryProfile[[3]] <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, sep=" ", stringsAsFactors=FALSE), silent=TRUE) # , quote=""
		queryProfile[[4]] <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, sep=";", stringsAsFactors=FALSE), silent=TRUE) # , quote=""
		err <- sapply(queryProfile, class)
		queryProfile <- queryProfile[err=="data.frame"]
		queryProfile <- queryProfile[sapply(queryProfile, dim)[2,] == max(sapply(queryProfile, dim)[2,])][[1]]
	    }
	    
# 	    queryProfile1 <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE))
# 	    queryProfile2 <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, quote="", sep=",", stringsAsFactors=FALSE))
# 	    queryProfile3 <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, quote="", sep=" ", stringsAsFactors=FALSE))
# 	    queryProfile4 <- try(read.table(file = paste0(userUploadedProfilePath,"/", exp), header=TRUE, quote="", sep=";", stringsAsFactors=FALSE))
# 	    if (dim(queryProfile2)[2] >= dim(queryProfile1)[2]) queryProfile <- queryProfile2 else queryProfile <- queryProfile1
# 	    if (dim(queryProfile3)[2] >= dim(queryProfile)[2]) queryProfile <- queryProfile3 else queryProfile <- queryProfile
# 	    if (dim(queryProfile4)[2] >= dim(queryProfile)[2]) queryProfile <- queryProfile4
	    
	    entries <- dim(queryProfile)[1]
	    if(entries == 0){
		    res$Remark <- "empty file"
		    return(res)
	    }
    # 			if(colnames(queryProfile)[1:4] == c("ID", "PROBE", "SYMBOL", "DESCRIPTION")) queryProfile <- queryProfile[, c("ID", "Diff_Exp", "Pvalues")]
# 	    if(identical(colnames(queryProfile)[1:4], c("ID_geneid", "PROBE", "Name_GeneSymbol", "DESCRIPTION"))) {
	    if(all(colnames(queryProfile)[1:4] == c("PROBE", "ID_geneid", "Name_GeneSymbol", "Description"))) {
		if(colnames(queryProfile)[(ncol(queryProfile)-1)] != "Value_LogDiffExp") {
# 			return(res) 
			queryProfile <- queryProfile[, c(2, ncol(queryProfile))]
		} else {
# 			queryProfile <- queryProfile[, c("ID", "SYMBOL", "Diff_Exp", "Pvalues")]
			queryProfile <- queryProfile[, c(2, 3, (ncol(queryProfile)-1), ncol(queryProfile))]
		}
	    }
	    
	    if(colnames(queryProfile)[1] == "signatureID") queryProfile <- queryProfile[,-1]
	    if("PROBE" %in% colnames(queryProfile)) queryProfile <- queryProfile[,-which(colnames(queryProfile)=="PROBE")]
	    if(dim(queryProfile)[2] == 4) {
		colnames(queryProfile) <- c("geneID","geneName","coefficients","Pvals")
		queryProfile <- queryProfile[!is.na(as.integer(queryProfile$geneID)), ]
		res$Remark <- "Done"
	    }
	    if(exp != "UpDn") res$notFound2 <- queryProfile[,1]
	    
# 	    if(dim(queryProfile)[2] == 4) queryProfile <- queryProfile[, c(1,3,4)]
	    if(dim(queryProfile)[2] == 3) colnames(queryProfile) <- c("geneID","coefficients","Pvals")
##	    if(dim(queryProfile)[2] == 3 & sum(is.na(as.integer(queryProfile[,2]))) > length(queryProfile[,2])/2) colnames(queryProfile) <- c("geneID","geneName","coefficients")
	    if(dim(queryProfile)[2] == 2) colnames(queryProfile) <- c("geneID","coefficients")
	    if(sum(is.na(as.integer(queryProfile[,1]))) > length(queryProfile[,1])/2) { ## sample with names only.  ## sum(sapply(as.integer(queryProfile[,1]), is.na)) != dim(queryProfile)[1]
		colnames(queryProfile)[1] <- "geneName"
		symid <- symbol2geneid(paste(queryProfile$geneName, collapse=","), org=org, debuging=debuging)
		symid <- symid[!duplicated(symid$Symbol), , drop=FALSE]
		si <- match(queryProfile$geneName, symid$Symbol)
		if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(geneID=symid$GeneID[si], geneName=queryProfile$geneName, coefficients=queryProfile$coefficients, Pvals=queryProfile$Pvals, stringsAsFactors=FALSE)
		if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(geneID=symid$GeneID[si], geneName=queryProfile$geneName, coefficients=queryProfile$coefficients, stringsAsFactors=FALSE)
# 		queryProfile <- as.data.frame(cbind(geneID=symid$GeneID[si], queryProfile), stringsAsFactors=FALSE)
		queryProfile <- queryProfile[!is.na(queryProfile$geneID), , drop=FALSE]
# 		if ("Pvals" %in% colnames(queryProfile)) {
# 		    queryProfile <- queryProfile[order(queryProfile$Pvals, decreasing=FALSE), , drop=FALSE]
# 			} else queryProfile <- queryProfile[order(queryProfile$geneID), , drop=FALSE]
		res$Remark <- "Done"
	    } else if(dim(queryProfile)[2] == 3 & sum(is.na(as.integer(queryProfile[,2]))) > length(queryProfile[,2])/2) { ## sample with names only.  ## sum(sapply(as.integer(queryProfile[,1]), is.na)) != dim(queryProfile)[1]
		colnames(queryProfile) <- c("geneID","geneName","coefficients")
# 		symid <- symbol2geneid(paste(queryProfile$geneName, collapse=","))
# 		symid <- symid[!duplicated(symid$Symbol), , drop=FALSE]
# 		si <- match(queryProfile$geneName, symid$Symbol)
# 		if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(geneID=symid$GeneID[si], geneName=queryProfile$geneName, coefficients=queryProfile$coefficients, Pvals=queryProfile$Pvals, stringsAsFactors=FALSE)
# 		if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(geneID=symid$GeneID[si], geneName=queryProfile$geneName, coefficients=queryProfile$coefficients, stringsAsFactors=FALSE)
# 		queryProfile <- as.data.frame(cbind(geneID=symid$GeneID[si], queryProfile), stringsAsFactors=FALSE)
		queryProfile <- queryProfile[!is.na(queryProfile$geneID), , drop=FALSE]
# 		if ("Pvals" %in% colnames(queryProfile)) {
# 		    queryProfile <- queryProfile[order(queryProfile$Pvals, decreasing=FALSE), , drop=FALSE]
# 			} else queryProfile <- queryProfile[order(queryProfile$geneID), , drop=FALSE]
		res$Remark <- "Done"
	    } else { 									## sample with IDs only.
		queryProfile <- queryProfile[!is.na(as.integer(queryProfile$geneID)), ,drop=FALSE]
		idsym <- geneid2symbol(paste(queryProfile$geneID, collapse=","), debuging=debuging)
		idsym <- idsym[!duplicated(idsym$GeneID), , drop=FALSE]
		si <- match(queryProfile$geneID, idsym$GeneID)
		if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(geneID=queryProfile$geneID, geneName=idsym$Symbol[si], coefficients=queryProfile$coefficients, Pvals=queryProfile$Pvals, stringsAsFactors=FALSE)
		if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(geneID=queryProfile$geneID, geneName=idsym$Symbol[si], coefficients=queryProfile$coefficients, stringsAsFactors=FALSE)
# 		queryProfile <- as.data.frame(cbind(geneID=queryProfile$geneID, geneName=idsym$Symbol[si], queryProfile[,-1, drop=FALSE]), stringsAsFactors=FALSE)
		queryProfile <- queryProfile[!is.na(queryProfile$geneID), , drop=FALSE]

		res$Remark <- "Done"
	    }
# # 	    if((dim(queryProfile)[2] %in% c(2,3)) & sum(is.na(as.integer(queryProfile[,1]))) > length(queryProfile[,1])/2) {  ## sum(sapply(as.integer(queryProfile[,1]), is.na)) != dim(queryProfile)[1]
# # 		colnames(queryProfile)[1] <- "geneName"
# # 		symid <- symbol2geneid(paste(queryProfile$geneName, collapse=","))
# # 		symid <- symid[!duplicated(symid$Symbol), , drop=FALSE]
# # 		si <- match(queryProfile$geneName, symid$Symbol)
# # 		queryProfile <- as.data.frame(cbind(geneID=symid$GeneID[si], queryProfile), stringsAsFactors=FALSE)
# # 		queryProfile <- queryProfile[!is.na(queryProfile$geneID), , drop=FALSE]
# # 		if ("Pvals" %in% colnames(queryProfile)) {
# # 		    queryProfile <- queryProfile[order(queryProfile$Pvals, decreasing=FALSE), , drop=FALSE]
# # 			} else queryProfile <- queryProfile[order(queryProfile$geneID), , drop=FALSE]
# # 		res$Remark <- "Done"
# # 	    }
# # 	    if (dim(queryProfile)[2] == 2) { #& sum(is.na(as.integer(queryProfile[,2]))) > length(queryProfile[,1])/2) {
# # 		queryProfile <- queryProfile[!is.na(as.integer(queryProfile$geneID)), ,drop=FALSE]
# # 		idsym <- geneid2symbol(paste(queryProfile$geneID, collapse=","))
# # 		idsym <- idsym[!duplicated(idsym$GeneID), , drop=FALSE]
# # 		si <- match(queryProfile$geneID, idsym$GeneID)
# # 		queryProfile <- as.data.frame(cbind(geneID=queryProfile$geneID, geneName=idsym$Symbol[si], queryProfile[,-1, drop=FALSE]), stringsAsFactors=FALSE)
# # 		queryProfile <- queryProfile[!is.na(queryProfile$geneID), , drop=FALSE]
# # 		queryProfile <- queryProfile[order(queryProfile$geneID), , drop=FALSE]
# # 		res$Remark <- "Done"
# # 	    }
# 	    if ((dim(queryProfile)[2] == 3) & sum(is.na(as.integer(queryProfile[,2]))) < length(queryProfile[,1])/2) {
# 		res$Remark <- "Done"
# 		}
	    
	    queryProfile <- unique(queryProfile)
	    if(dim(queryProfile)[2] != 2 & dim(queryProfile)[2] != 3 & dim(queryProfile)[2] != 4){
		    res$Remark <- "wrong format"
		    return(res)
	    }
	    
# 	    queryProfile <- na.omit(queryProfile)

## new addition to get non human number of genes
	    homoloArr <- getHomologousGenes(as.character(queryProfile[,1]), org=org)
# # # 	    shorthomoloArr <- getHomologousGenes(as.character(queryProfile[1:min(100, nrow(queryProfile)),1]), org=org)
	    
	    if(length(!is.na(homoloArr[,3])) == 0){
		    res$Remark <- paste0("No homolo genes found from ", org)
		    return(res)
	    }
	    queryProfile[,1:2] <- homoloArr[match(queryProfile[,1], homoloArr[,1]),3:4]
	    queryProfile <- queryProfile[which(!is.na(queryProfile[,1])),]
# 	    rownames(queryProfile) <- queryProfile$geneID
## end new addition
	    
	    if(exp!="UpDn") {
		if ("Pvals" %in% colnames(queryProfile)) {
		    queryProfile <- queryProfile[order(queryProfile$Pvals, decreasing=FALSE), , drop=FALSE]
		    if (!is.null(topGenes)) queryProfile <- queryProfile[1:topGenes, , drop=FALSE]
		} else {
		    queryProfile <- queryProfile[order(abs(queryProfile$coefficients), decreasing=TRUE), , drop=FALSE]
		    if (!is.null(topGenes)) queryProfile <- queryProfile[1:topGenes, , drop=FALSE]
		}
	    }
	    
	    #####Generating sessionID
	    # Note: in other cases, sessionID is written by plot_with_properties_eset function in mod_Rserve_prod.R - Nick
	    sessionID <- paste(c(strsplit(date(),split=" ")[[1]],as.integer(runif(1)*10e6)),sep="",collapse="_")
	    sessionID <- gsub(":","_", sessionID)
	    sessionID <- gsub("__","_", sessionID)
	    data(LibGenes)
	    commonGenes <- as.data.frame(lapply(LibGenes, function(j) length(intersect(j, homoloArr[,3])))) # queryProfile$geneID
##	    shortqueryProfile <- queryProfile[1:min(100, nrow(queryProfile)),]
# # # 	    shortcommonGenes <- as.data.frame(lapply(LibGenes, function(j) length(intersect(j, shorthomoloArr[,2])))) # shortqueryProfile$geneID
	    
	    if (write) {
		if(all(abs(queryProfile[,3])==1)) {
		    shortqueryProfile <- queryProfile
		} else {
		    shortqueryProfile <- queryProfile[1:min(100, nrow(queryProfile)),]
		}
		shortfileName <- paste0("processedSig_short_", sessionID,".xls")
		fileName <- paste0("processedSig_", sessionID,".xls")
		if (dim(queryProfile)[2] == 3) {
		    colnames(queryProfile) <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
		    colnames(shortqueryProfile) <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
		} else {
		    colnames(queryProfile) <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
		    colnames(shortqueryProfile) <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
		}
		write.table(queryProfile, file=paste0(path_to_write, fileName), quote=TRUE, sep="\t", row.names=FALSE)
		write.table(shortqueryProfile, file=paste0(path_to_write, shortfileName), quote=TRUE, sep="\t", row.names=FALSE)
		res$fileName <- fileName
		res$shortfileName <- shortfileName
	    }
	    
	    res$sessionID <- sessionID
# 	    res$Remark <- "Done"
	    res$entries <- entries
	    
	    geneinfo <- get(data(geneinfo))
	    ids <- geneinfo[geneinfo$GeneID %in% res$notFound2, ]
	    syms <- geneinfo[geneinfo$Symbol %in% res$notFound2, ]
	    id <- F
	    if(dim(ids)[1] > dim(syms)[1]) {
		tmp <- ids 
		id <- T
	    } else tmp <- syms
	    if(id) {
		##homo <- getHomologousGenes(as.character(tmp$GeneID), org=org)
		res$notFound1 <- setdiff(res$notFound2, tmp$GeneID)
		##res$notFound2 <- tmp$GeneID[tmp$TaxID != org2tax(org)]
		res$notFound2 <- homoloArr[is.na(homoloArr$homoloGeneID), "geneID"]
	    } else {
		##homo <- getHomologousGenes(as.character(tmp$GeneID), org=org)
		res$notFound1 <- setdiff(res$notFound2, tmp$Symbol)
		##es$notFound2 <- tmp$GeneID[tmp$TaxID != org2tax(org)]
		res$notFound2 <- homoloArr[is.na(homoloArr$homoloGeneID), "Symbol"]
# 		res$notFound2 <- tmp$GeneID[!match(tmp$GeneID, homo$geneID)]
	    }

	    
	    if(sum(is.na(as.integer(res$notFound))) > dim(queryProfile)[1]/2) {
		res$notFound <- setdiff(res$notFound, queryProfile$Name_GeneSymbol)
	    } else {
		res$notFound <- setdiff(res$notFound, queryProfile$ID_geneid)
	    }
# # # 	    res$shortSig <- shortqueryProfile
	    res$Sig <- queryProfile
	    res$commonGenes <- commonGenes
# # # 	    res$shortcommonGenes <- shortcommonGenes
	    res$foundGenes <- dim(queryProfile)[1]
	    res$foundUniqueGenes <- length(unique(queryProfile[,1]))
	    if (ncol(queryProfile)==3) res$Pval <- FALSE
	    res
#@!    }
}


## How to get all gene IDs
# tmp <- list(LIB_1=NA, LIB_2=NA, LIB_3=NA, LIB_5=NA, LIB_6=NA, LIB_7=NA, LIB_8=NA, LIB_10=NA)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/gdsDiffExp.RData")
# tmp$LIB_1 <- rownames(gdsDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/cmapDiffExp.RData")
# tmp$LIB_2 <- rownames(cmapDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/encodeScores.RData")
# tmp$LIB_3 <- rownames(encodeScores)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/lincscpDiffExp.RData")
# tmp$LIB_5 <- rownames(lincscpDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/lincskdDiffExp.RData")
# tmp$LIB_6 <- rownames(lincskdDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/lincsshDiffExp.RData")
# tmp$LIB_7 <- rownames(lincsshDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/lincsRNAseqDiffExp.RData")
# tmp$LIB_8 <- rownames(lincsRNAseqDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/ctrsDiffExp.RData")
# tmp$LIB_10 <- rownames(ctrsDiffExp)
# 
# ## noticed new 2018 generated signature files have different rownames (gene ids) order!!!
# load("/opt/raid10/lincs/ilincsDatabase/signatures/2018/lincsoeDiffExp.RData")
# tmp$LIB_11 <- rownames(lincsoeDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/2018/drugmatrixDiffExp.RData")
# tmp$LIB_12 <- rownames(drugmatrixDiffExp)
# 
# load("/opt/raid10/lincs/ilincsDatabase/signatures/2018/lincskdDiffExp_bing.RData")
# tmp$BING <- rownames(lincskdDiffExp) ## this is a different file than above!
# 
# LibGenes <- tmp
# save(LibGenes, file="/opt/raid10/genomics/mehdi/ilincs/gitHub/ilincsR/data/LibGenes.rda", compress=FALSE)


