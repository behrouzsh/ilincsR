#' A function for visualizing, analysing and further investigating on a next-gene experiment.
#'
#' This function allows you to run multiple analysis and is one of the main work flows in GenomicsPortals.
#' @param exp This argument is the experiment name, defined as a string corresponding to GenomicsPortals experiment name. The default is set to NULL.
#' @param path_to_write This necessary parameter specifies the path which user wants to save the
#'	output files (heatmaps, TreeView files and some other output files) to be saved in.
#' @param libName libName defines which precomputed signature library to calculate concordances against.
#' @param userUploadedProfile Is required if user wants to use his/her own provided gene Profiles in the analysis. It should be in the form of a string containig the path to Profile text file.
#'	This text file should contain geneIDs in the first column, Coefficients in the second and pValues in the third column (It's also possible to submit the file without pValues).
#' @param queryGenelist This is a string parameter consisting of user provided gene IDs separated by commas to be used in querying the dataset or experiment.
#' @param write write argument is used (in downstreem function "parse_res") to save the SignatureTable for further uses.
#'	 It should be set to TRUE, as it is by default, to save the SignatureTable in the form of txt.
#' @param homoloGenes homoloGenes To find human homologus genes if signature is not from human.
#' @param cutoff Is the p-value cutoff to filter out concordance results.
#' @param debuging for debugging purposes only!.
#' @param sigDB Which precomputed signature database to use, before or after 2018.
#' @param algorithm Different data source and connectivity methods including c("Cor_all_exp", "Cor_all_plog", "Cor_updn_exp", "Cor_updn_plog", "wtCor_all_exp").
#' @param chunk Used in splitting sig data for ENCODE only.
#' @param ncors Used in multi thread calculation of ENCODE sig data only.
#' @param org Organism used for mapping genes.
#' @param authors M. Fazel-Najafabadi
#' @keywords exp, userUploadedProfile, algorithm
#' @export
#' @examples
#' ## not run
#' # For uploading a signature and finding concordances to a library:
#' res <- findConcordances(exp="sample.tsv", path_to_write="/patho/to/your/output/directory/", libName="LIB_5",
#'	userUploadedProfile="/patho/to/your/output/directory/")
#'
#' ## not run
#' @export

findConcordancesBingX <- function(exp=NULL, path_to_write, libName=NULL, userUploadedProfile=NULL, queryGenelist=NULL, write=FALSE,
				   homoloGenes=TRUE, cutoff=0.05, logPcut=10, debuging=FALSE, sigDB="ilincs_sigs", algorithm="wtCor_all_exp", topSigs=FALSE, metadata=FALSE, chunk=500, ncors=10, org="Hs")
{

  time0 <- Sys.time()

## Mehdi: log

#        print(paste("exp",exp))
#        print(paste("path_to_write",path_to_write))
#        print(paste("libName",libName))
#        print(paste("userUploadedProfile",userUploadedProfile))
#        print(paste("queryGenelist",queryGenelist))
#	 print(paste("before get pdata"))

	print(paste0("libName is ", libName))
	if(libName %in% c("LIB_2","LIB_3","LIB_8","LIB_9")) topSigs <- FALSE

	res <- list(sessionID=NA, gpgene=NA, gpprobe=NA, Remark=NA, NoOfGenes=NA, NoOfProbes=NA, concordanceTable=NA)
	servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
	if(servSet$sigDB == "ilincs_new") algorithm <- "Cor_all_exp"

	## Mehdi: adding signatures meta data
        if(metadata) {
		metacols <- switch(libName,
			    LIB_1="datasetID,factor,Level1,Level2",
			    LIB_2="compound,lincsPertID,concentration,cellLine",
			    LIB_3="treatment,antibodyTarget,cellLine",
			    LIB_5="compound,lincsPertID,integratedGeneTargets,concentration,cellLine,time",
			    LIB_6="treatment,perturbagenID,cellLine,time",
			    LIB_7="treatment,perturbagenID,cellLine,time",
			    LIB_8="datasetID,factor,Level1,Level2",
			    LIB_9="compound,lincsPertID,concentration,cellLine,time",
			    LIB_10="perturbagenID,compound,lincsPertID,cellLine",
			    LIB_11="treatment,perturbagenID,cellLine,time",
			    LIB_12="tissue,compound,concentration,organism", #,lincsPertID
			    LIB_13="datasetID,factor,Level1,Level2",
			    LIB_14="compound,integratedGeneTargets,concentration,cellLine,time"
		)
	} else {
		metacols <- switch(libName,
			    LIB_1="signatureID",
			    LIB_2="signatureID",
			    LIB_3="signatureID",
			    LIB_5="lincsSigID,signatureID",
			    LIB_6="lincsSigID,signatureID",
			    LIB_7="lincsSigID,signatureID",
			    LIB_8="signatureID",
			    LIB_9="signatureID",
			    LIB_10="signatureID",
			    LIB_11="lincsSigID,signatureID",
			    LIB_12="signatureID",
			    LIB_13="signatureID",
			    LIB_14="signatureID"
		)
	}

        if(metadata) {
		cnames <- switch(libName,
			    LIB_1=c("signatureID", "similarity", "pValue", "nGenes", "datasetID","factor","Level1","Level2"),
			    LIB_2=c("signatureID", "similarity", "pValue", "nGenes", "compound", "lincsPertID", "concentration","cellLine"),
			    LIB_3=c("signatureID", "similarity", "pValue", "nGenes", "treatment","antibodyTarget","cellLine"),
			    LIB_5=c("signatureID", "similarity", "pValue", "nGenes", "compound", "lincsPertID", "GeneTargets", "concentration", "cellLine", "time"),
			    LIB_6=c("signatureID", "similarity", "pValue", "nGenes", "treatment","perturbagenID","cellLine","time"),
			    LIB_7=c("signatureID", "similarity", "pValue", "nGenes", "treatment","perturbagenID","cellLine","time"),
			    LIB_8=c("signatureID", "similarity", "pValue", "nGenes", "datasetID","factor","Level1","Level2"),
			    LIB_9=c("signatureID", "similarity", "pValue", "nGenes", "compound", "lincsPertID", "concentration", "cellLine", "time"),
			    LIB_10=c("signatureID", "similarity", "pValue", "nGenes", "perturbagenID","compound", "lincsPertID", "cellLine"),
			    LIB_11=c("signatureID", "similarity", "pValue", "nGenes", "treatment","perturbagenID","cellLine","time"),
			    LIB_12=c("signatureID", "similarity", "pValue", "nGenes", "tissue", "compound", "concentration","organism"), #, "lincsPertID"
			    LIB_13=c("signatureID", "similarity", "pValue", "nGenes", "datasetID","factor","Level1","Level2"),
			    LIB_14=c("signatureID", "similarity", "pValue", "nGenes", "compound", "GeneTargets", "concentration", "cellLine", "time")
		)
	} else {
		cnames <- switch(libName,
			    LIB_1=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_2=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_3=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_5=c("lincsSigID", "signatureID", "similarity", "pValue", "nGenes"),
			    LIB_6=c("lincsSigID", "signatureID", "similarity", "pValue", "nGenes"),
			    LIB_7=c("lincsSigID", "signatureID", "similarity", "pValue", "nGenes"),
			    LIB_8=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_9=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_10=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_11=c("lincsSigID", "signatureID", "similarity", "pValue", "nGenes"),
			    LIB_12=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_13=c("signatureID", "similarity", "pValue", "nGenes"),
			    LIB_14=c("signatureID", "similarity", "pValue", "nGenes")
		)
	}
	if(exp=="geneList") {
	    queryGenelist <- userUploadedProfile[[1]] # unlist(strsplit(queryGenelist,","))
	    usedGenes <- length(unique(queryGenelist))
	    enrichTable <- computeGenelistEnrichment(queryGenelist=queryGenelist, libName=libName, debuging=debuging)
	    if(dim(na.omit(enrichTable))[1] > 2){ remark <- "Done"
	    } else{ remark <- "Error in Random Set or No enriched genes found"
		enrichTable <- na.omit(enrichTable)}
	    sessionID <- generateSID()
	    enrichTablePath <- paste(path_to_write,"enrichTable_",sessionID,".xls",sep="")
	    sigScores <- data.frame(signature=enrichTable[, 1], score=round(enrichTable[, 3],4), stringsAsFactors=FALSE)
	    write.table(enrichTable[,-2],file=paste(path_to_write,"enrichTable_",sessionID,".xls",sep=""), sep="\t", quote=TRUE, row.names=F, col.names=F, append=T)
	    sigmeta <- getSignatureMeta(prop=metacols, signatures=paste(sigScores$signature, collapse=","), debuging=debuging, sigDB=servSet$sigDB)
	    sigmatch <- match(sigScores$signature, sigmeta$signatureID)
	    sigmeta <- sigmeta[sigmatch,]
	    sigmeta <- sigmeta[,-1]
	    if (libName=="LIB_1" | libName=="LIB_8" | libName=="LIB_13") {
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:4]) #else sigScores <- cbind(sigScores, sigmeta[,"lincsSigID"])
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "datasetID","factor","Level1","Level2") else colnames(sigScores) <- c("signatureID", "zScores")
	    } else if  (libName=="LIB_2") {
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:3]) #else sigScores <- cbind(sigScores, sigmeta[,"lincsSigID"])
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "compound","concentration","cellLine") else colnames(sigScores) <- c("signatureID", "zScores")
	    } else if  (libName=="LIB_12") {
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:4]) #else sigScores <- cbind(sigScores, sigmeta[,"lincsSigID"])
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "tissue", "compound","concentration","organism") else colnames(sigScores) <- c("signatureID", "zScores")
	    } else if (libName=="LIB_3") {
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:3]) #else sigScores <- cbind(sigScores, sigmeta[,"lincsSigID"])
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "treatment","antibodyTarget","cellLine") else colnames(sigScores) <- c("signatureID", "zScores")
	    } else if (libName=="LIB_5") {
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:6]) else sigScores <- cbind(sigmeta[,"lincsSigID"], sigScores)
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "compound","perturbagenID","ProteinTargets","concentration","cellLine","time") else colnames(sigScores) <- c("lincsSigID", "signatureID", "zScores")
	    } else if (libName=="LIB_10") {
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:3]) #else sigScores <- cbind(sigScores, sigmeta[,"lincsSigID"])
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "perturbagenID", "compound", "cellLine") else colnames(sigScores) <- c("signatureID", "zScores")
	    } else {#	if (libName="LIB_6" | libName="LIB_7")
		if(metadata) sigScores <- cbind(sigScores, sigmeta[,1:4]) else sigScores <- cbind(sigmeta[,"lincsSigID"], sigScores)
		if(metadata) colnames(sigScores) <- c("tableNames", "zScores", "treatment","perturbagenID","cellLine","time") else colnames(sigScores) <- c("lincsSigID", "signatureID", "zScores")
	    }
	    res <- list(usedGenes=usedGenes , sessionID=sessionID, Remark=remark, enrichTablePath=enrichTablePath, sigScores=sigScores)
	    return(res)

	} else if(exp=="UpDn") { sessionID <- generateSID()
	    res$Remark <- "Done"
	    res$sessionID <- sessionID
	    upgenes <- userUploadedProfile[[1]] #$genesUp #
	    dngenes <- userUploadedProfile[[2]] #$genesDown #
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
	} else if(exp=="fullSig") {
	    sessionID <- generateSID()
	    res$Remark <- "Done"
	    res$sessionID <- sessionID
	    queryProfile <- userUploadedProfile[[1]]
	} else if(!is.null(userUploadedProfile)) {
		if (length(grep("processedSig_", exp)) == 1) {
		  queryProfile0 <- try(read.table(file = paste(userUploadedProfile,"/", exp, sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)) # , quote=""
		} else {
		  queryProfile0 <- uploadedSigProcess(exp, userUploadedProfile, write=write, debuging=debuging)$Sig
		}
		if (dim(queryProfile0)[2] == 3) {
		    colnames(queryProfile0) <- c("geneID","geneName","coefficients")
		} else {
		    colnames(queryProfile0) <- c("geneID","geneName","coefficients","Pvals")
		}
		queryProfile <- queryProfile0[, -2]
		sessionID <- generateSID()
		res$Remark <- "Done"
		res$sessionID <- sessionID
	} else {
		res$Remark <- "No profile uploaded."
		return(res)
	}




# 	  homoloArr <- getHomologousGenes(as.character(queryProfile[,1]), org=org)
# 	  if(length(!is.na(homoloArr[,3])) == 0){
# 		  res <- getDefaultErrorDataframe("No homolo genes")
# 	  }
# 	  queryProfile[,1] <- homoloArr[match(queryProfile[,1], homoloArr[,1]),3]
# 	queryProfile <- queryProfile[which(!is.na(queryProfile[,1])),]
# 	if(length(unique(queryProfile[,"geneID"])) < dim(queryProfile)[1]) {
# 	    if(dim(queryProfile)[2] == 3) unqIdPvals <- split(queryProfile[,"Pvals"], queryProfile[,"geneID"])
# 	    if(dim(queryProfile)[2] == 3) mergePvals <- unlist(lapply(unqIdPvals, function(x) exp(mean(log(x), na.rm = TRUE)))) ## Mehdi: I changed it so we dont loose genes with p-val=NA
# 	    unqIdCoeffs <- split(queryProfile[,"coefficients"], queryProfile[,"geneID"])
# 	    mergeCoeffs <- unlist(lapply(unqIdCoeffs, function(x) mean(as.numeric(x), na.rm = TRUE)))
# 	    if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(geneID=names(unqIdCoeffs), coefficients=mergeCoeffs, Pvals=mergePvals, stringsAsFactors=FALSE)
# 	    if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(geneID=names(unqIdCoeffs), coefficients=mergeCoeffs, Pvals=1, stringsAsFactors=FALSE)
# 	}
# 	if(dim(queryProfile)[2] == 2) queryProfile$Pvals <- 1
# 	if(all(abs(queryProfile$coefficients) == 1)) algorithm <- "Cor_all_exp"
# 	rownames(queryProfile) <- queryProfile$geneID ## Mehdi: I added this later on
#         mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=servSet$sigDB, host=servSet$host, port = servSet$port, password="public")
#         sql <- paste0("select diffexpTablePath,pvaluesTablePath,scoresTablePath,probsTablePath from signatureLibraries where libraryID = '", libName, "'")
#         dt <- DBI::dbGetQuery(mycon, sql)
#         libraryPath <- dt[1,"diffexpTablePath"]
#         libraryPvalues <- dt[1,"pvaluesTablePath"]
# 	libraryScore <- dt[1,"scoresTablePath"]
# 	libraryProb <- dt[1,"probsTablePath"]
#         DBI::dbDisconnect(mycon)
# if (libName=="LIB_3"){
#     refProfiles <- get(load(libraryScore)); refProbs <- get(load(libraryProb))
#
# } else {
#
#     if(servSet$sigDB=="ilincs_new") {
# 	refProfiles <- get(load(libraryPath))
#     } else {
# 	if(algorithm=="Cor_all_exp" & servSet$sigDB != "ilincs_new") {
# 	    refProfiles <- get(load(libraryPath))
# 	}
# 	if(algorithm %in% c("Cor_all_plog", "Cor_updn_plog")) {
# 	    libraryPvalues <- sub("2018/", "2018/signed_pLogs/", libraryPvalues)
# 	    libraryPvalues <- sub("PValues", "Plogs", libraryPvalues)
# 	    refProfiles <- get(load(libraryPvalues))
# 	}
# 	if(algorithm %in% c("Cor_updn_exp", "wtCor_all_exp")) {
# 	    refProfiles <- get(load(libraryPath))
# 	    libraryPvalues <- sub("2018/", "2018/signed_pLogs/", libraryPvalues)
# 	    libraryPvalues <- sub("PValues", "Plogs", libraryPvalues)
# 	    refProbs <- get(load(libraryPvalues))
# 	}
#
#     }
# }
# 	print(paste(Sys.time() - time0, ": LIB path: ", libraryPath))
#
# 	refgenes <- get(data(LibGenes))[[libName]]
# 	querygenes <- as.character(queryProfile[,1]) ## Mehdi: why it was integer?!!!
# 	genes <- intersect(refgenes, querygenes)
# 	queryData <- queryProfile[genes, ]
#
#
# 	if (libName=="LIB_3" & dim(queryProfile)[2] == 3){
# 		refData <- refProfiles[genes, ,drop=FALSE]
# 		no <- apply(refData,2,function(i) all(is.na(i)))
# 		refData <- refData[, !no ,drop=FALSE]
# 		refProbsData <- refProbs[genes, ,drop=FALSE]
# 		refProbsData <- refProbsData[, !no ,drop=FALSE]
# 	} else {
# 		refData <- refProfiles[genes, ,drop=FALSE]
# 		no <- apply(refData,2,function(i) all(is.na(i)))
# 		refData <- refData[, !no ,drop=FALSE]
# 		if(algorithm == "Cor_updn_exp") {
# 		    u <- d <- ifelse(nrow(refData) > 2500, 500, 100)
# 		    refProbsData <- refProbs[genes, !no,drop=FALSE]
# 		    refData <- setUpDn(refData, refProbsData, u, d, 0)
# 		}
# 		if(algorithm == "Cor_updn_plog") {
# 		    u <- d <- ifelse(nrow(refData) > 2500, 500, 100)
# 		    refData <- setUpDn(refData, refData, u, d, 0)
# 		}
# 		if(algorithm=="wtCor_all_exp") {
# 		    refProbsData <- refProbs[genes, !no,drop=FALSE]
# 		}
# 	}
# 	if (dim(queryData)[1] == 0) {
# 		res$Remark <- "Filtered zero genes" # no overlapping genes between query and reference
# 		print("********************* Filtered zero genes!")
# 		return(res)
# 	}
# 	if (dim(queryData)[1] == 1) {
# 		res$Remark <- paste0("Filtered one. Only geneID ", rownames(queryData), " is in common.") # no overlapping genes between query and reference
# 		print(paste0("********************* Filtered one. Only geneID ", rownames(queryData), " is in common."))
# 		return(res)
# 	}
# 	if (sd(queryData[,2]) == 0) {
# 		res$Remark <- "Uniform expression for common genes in signature. No correlation!" # no overlapping genes between query and reference
# 		print("********************* Uniform expression in signature for common genes. No correlation!")
# 		return(res)
# 	}
# 	if(libName=="LIB_3"){
# 	    if(dim(queryProfile)[2] == 3) {
# 		queryDataDown <- ifelse(queryData[,2]<0,queryData[,3]/2, 1-queryData[,3]/2)
# 		queryDataDown[queryDataDown<1e-10] <- 1e-10
# 		queryDataDown <- data.frame(queryData[,1],queryDataDown)
# 		queryDataUp <- ifelse(queryData[,2]>0,queryData[,3]/2, 1-queryData[,3]/2)
# 		queryDataUp[queryDataUp<1e-10] <- 1e-10
# 		queryDataUp <- data.frame(queryData[,1],queryDataUp)
# 		 grsZAll <- parallel::mclapply(1:dim(refData)[2],function(j) {
# 				print(j)
#  				combineData <- na.omit(data.frame(EntrezID=queryDataUp[,1],Score=refData[,j],Prob=refProbsData[,j],queryDataUp[,2],queryDataDown[,2],stringsAsFactors=F))
# 				grsZUp <- batchTregGRS(query.p=combineData[,c(1,4)],reference.p=combineData[,c(1,2,3)] , na.rm = TRUE, estimateNullDistr = F, nullDistrQuantiles = c(0.9, 0.95, 0.99),
#                                		nullDistrN = 100, tolerateWarnings=TRUE, pAdjust.method.query=NULL, pAdjust.method.reference=NULL, lambda=0.005,plotRescaled=F)[1:2]
#      				grsZDown <- batchTregGRS(query.p=combineData[,c(1,5)],reference.p=combineData[,c(1,2,3)] , na.rm = TRUE, estimateNullDistr = F, nullDistrQuantiles = c(0.9, 0.95, 0.99),
#                                		nullDistrN = 100, tolerateWarnings=TRUE, pAdjust.method.query=NULL, pAdjust.method.reference=NULL, lambda=0.005,plotRescaled=F)[1:2]
#   #                 		grsZ <- ifelse( abs(grsZUp)> abs(grsZDown),grsZUp, grsZDown*(-1))
#                    		grsP <- ifelse( abs(grsZUp[[1]])> abs(grsZDown[[1]]), grsZUp[[1]], ifelse(grsZDown[[1]] < 0, grsZDown[[1]], grsZDown[[1]]*(-1)))
#                    		grsZ <- ifelse( abs(grsZUp[[2]])> abs(grsZDown[[2]]), grsZUp[[2]], ifelse(grsZDown[[2]] < 0, grsZDown[[2]], grsZDown[[2]]*(-1)))
#    				list(grsZ, grsP)
#  			}, mc.cores=ncors)
#  		grsZAll <- as.data.frame(data.table::rbindlist(grsZAll))
# 		rownames(grsZAll) <- colnames(refData)
# 		ord <- order(abs(as.numeric(grsZAll[,1])),decreasing=TRUE)
# 		#ord <- order(unlist(grsZAll),decreasing=TRUE)
#                 grsZAll <- grsZAll[ord,]
# 		#return All significant or top 10
# 		sigComp <- which(abs(as.numeric(grsZAll[,1])) >= 6)
# 		if(length(sigComp) == 0) { sigComp <- c(1:min(10, nrow(grsZAll))) }
# 		nGenes <- t(!is.na(queryData[,1])) %*% (!is.na(refData))
# 		corTable <- data.frame(dataset=rownames(grsZAll)[sigComp],CC=grsZAll[sigComp,1], PV=grsZAll[sigComp,2], nGenes=nGenes[1,sigComp], stringsAsFactors=F)
# 	    } else{
# 		res$Remark <- "No P-values provided"
# 		return(res)
# 	    }
#
# 	} else {
# 		if(algorithm != "wtCor_all_exp") {
# 		    if(dim(queryData)[1] == 2) {
# 			tmpcor <- cor(queryData[,2], refData, use="pairwise.complete.obs")
# 			AllCor <- data.frame(dataset=colnames(tmpcor), CC=tmpcor[1,], PV=NA, nGenes=2, stringsAsFactors=FALSE)
# 		    } else {
# 			if(algorithm == "Cor_updn_exp") {
# 			    u <- d <- ifelse(nrow(queryData) > 2500, 500, 100)
# 			    queryData[,3] <- sign(queryData[,2]) * (-log10(queryData[, 3]))
# 			    queryData[,2] <- setUpDn(queryData[,2,drop=F], queryData[,3,drop=F], u, d, 0)
# 			}
# 			if(algorithm == "Cor_updn_plog") {
# 			    u <- d <- ifelse(nrow(queryData) > 2500, 500, 100)
# 			    queryData[,3] <- sign(queryData[,2]) * (-log10(queryData[, 3]))
# 			    queryData[,2] <- setUpDn(queryData[,3,drop=F], queryData[,3,drop=F], u, d, 0)
# 			}
# 			AllCor <- cor_with_pval(queryData[,2], refData, cortype="pearson")#use="pairwise.complete.obs")
# 		    }
# 		} else {

	library(iLincsCor)
	libs<-get("libs", envir=.pkgglobalenv)
  	x_lib <- libs[[libName]]
  	if (is.null(x_lib)) {
		return(NULL)
	}
	vec <- x_lib$read_input(paste(userUploadedProfile,"/", exp, sep=""))

            # queryData[is.na(queryData[,3]),3] <- 1
            # log10p <- - log10(queryData[,3])
            # log10p[log10p > logPcut] <- logPcut
            # wt <- abs(refProbsData) + log10p

	NoOfGenes <- vec$found_geneids

	# save(vec, file="mk_vec.RDS")

	AllCorVec <-x_lib$cor(vec,4)

	# save(AllCorVec, file="mk_AllCorVec.RDS")

	include <- which(abs(AllCorVec)>=0.2)
	AllCor<-data.frame(dataset=x_lib$signature_names[include], CC=AllCorVec[include], PV=0, nGenes=0)


		#   print(paste(Sys.time() - time0,"R START -- mehdi_wtcor"))
		#   x<-system.time({
		# 	queryData[is.na(queryData[,3]),3] <- 1
		# 	log10p <- - log10(queryData[,3])
		# 	log10p[log10p > logPcut] <- logPcut
		# 	wt <- abs(refProbsData) + log10p
		# 	AllCor <- mehdi_wtcor(queryData[,2], refData, wt, topSigs=topSigs)
		# 	AllCor <- data.frame(signatureID=colnames(AllCor), t(AllCor)[,c(1,4,3)], stringsAsFactors=F)
		#	  AllCor[,4] <- AllCor[,4] + 2
		#   })
		#   print(x)
		# 	print(paste(Sys.time() - time0,"R END -- mehdi_wtcor"))
		# }
		#colnames(AllCor) <- c("dataset", "CC", "PV", "nGenes")
		#ord <- order(abs(AllCor$CC), decreasing=TRUE)
		#AllCor <- AllCor[ord,]
		#AllCor <- AllCor[!is.na(AllCor$PV),,drop=F]
		# sigComp <- AllCor$PV <= cutoff & abs(AllCor$CC) >= 0.2
	  if (nrow(AllCor) == 0) {
	    ord <- order(abs(AllCorVec), decreasing=TRUE)
	    Top10CC <- head(AllCorVec[ord],n=10)
	    Top10names <- head(x_lib$signature_names[ord],n=10)

	    AllCor<-data.frame(dataset=Top10names, CC=Top10CC, PV=0, nGenes=0)
	  }
		# if(sum(sigComp, na.rm=T) == 0){ # no significant comparisons then return top 10
		# 	corTable <- AllCor[1:min(10,nrow(AllCor)),] # sigComp <- rep(T,10)
		# 	corTable <- corTable[!is.na(corTable$dataset),]
		# } else {
		# 	corTable <- AllCor[sigComp,]
		# 	corTable <- corTable[!is.na(corTable$dataset),]
		# 	corTable <- corTable[1:min((nrow(AllCor)%/%10),nrow(corTable)),]
		# }
	corTable <- AllCor
	# }

    if (metadata) {
	   sigmeta <- getSignatureMeta(signatures=paste(corTable$dataset, collapse=","), prop=metacols, debuging=debuging, sigDB=servSet$sigDB)
	   sigmatch <- match(corTable$dataset, sigmeta$signatureID)
	   sigmeta <- sigmeta[sigmatch,,drop=FALSE]
	   sigmeta <- sigmeta[,-1,drop=FALSE]
	   if(metadata) {
	       corTable <- cbind(corTable, sigmeta)
	   } else {
	       corTable <- cbind(sigmeta, corTable[,-1]) #[,c(5,2,3,4)]
	   }
	   colnames(corTable) <- cnames
    }


	res$concordanceTable <- corTable
	write.table(corTable,file=paste(path_to_write,"CorTable_",as.character(res[[1]]),".xls",sep=""), quote=TRUE, sep="\t", row.names=F, col.names=T, append=T)
        res[[(length(res)+1)]] <- paste(path_to_write,"CorTable_",as.character(res[[1]]),".xls",sep="")
        names(res)[length(res)] <- "corTablePath"
        res$NoOfGenes <- NoOfGenes # length(genes)
     	return(res)
}

