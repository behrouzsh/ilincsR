
#' A function for analysing and visualizing signatures.
#'
#' This function allows you to analyse and visualize precomputed signatures available in iLincs portal.
#' @param idList A list of signature names as a single string/character argument separated by commas.
#' @param libName not used so far!
#' @param path_to_write This parameter specifies the path which user wants to save the results.
#' @param queryID not used so far!
#' @param NoOfGenesForHeatmap The number of top genes per signature based on their significance level user wants to keep in the analysis.
#' @param glist A list of Entrez gene IDs as a single string/character argument 
#'	separated by commas user wants to investigate in the selected signatures.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param sigDB The signatures version available on iLincs portal. It can have values "new" or "old". The default is  "new".
#' @param author M. Fazel-Najafabadi
#' @keywords ...
#' @export 
#' @examples
#'	## not run
#'	idList <- "LINCSKD_43588,LINCSKD_100,LINCSKD_49086,LINCSCP_101"
#'	glist <- "1514,3553,8835,727,23635,23463,8870,2817,2195,3280,23200,10868,6774,2887,5525,3486,23588,
#'		3122,10513,3693,9019,5747,10444,22926,25976,9961,1029,9053,26511,593,1891,11031,28969,148022,
#'		6850,79174,23,5321,9842,64764,2263,7016,10123,10775,6118,622,80025,10058,10153,9128,6182,2553,
#'		9801,9805,10953,9217,4927,5257,1070,10318,3988,5696,90861,2958,23210,891,2050,84159,1282,1845,
#'		7750,1017"
#'	mastergroup <- master_GroupLincsAnalysis(idList=idList, path_to_write="your/path/here", glist=glist, sigDB="ilincs_sigs")
#'	# It can be easily noticed that all these signatures are present in the old version of signatures
#'	# but not in the new version.
#'	## end not run

proteinTargets <- function(pertList, path_to_write, debuging=FALSE, sigDB="ilincs_sigs") {
    
    sessionID <- generateSID()
    esetNameP <- paste0(sessionID, "_ProteinTargets")
    targets <- drugTargets(drugID = pertList, debuging = debuging, sigDB = sigDB)
    allTargets <- paste(targets$targets, collapse = ", ")
    allTargets <- unique(unlist(strsplit(allTargets, ", ")))
    geneInfo <- symbol2geneid(paste(allTargets, collapse = ","), debuging = debuging)
    geneInfo <- geneid2symbol(paste(geneInfo$GeneID, collapse = ","), description = T, debuging = debuging)
    rownames(geneInfo) <- geneInfo$Symbol
    
    
    targetsMatrix <- matrix(0, nrow = nrow(targets), ncol = nrow(geneInfo))
    rownames(targetsMatrix) <- targets$drugID
    colnames(targetsMatrix) <- geneInfo$Symbol
    
    for(i in rownames(targetsMatrix)) {
	targetsMatrix[i, unlist(strsplit(targets[targets$drugID == i, "targets"], ", "))] <- 10
    }
    
    met <- data.frame(drugID=rownames(targetsMatrix), stringsAsFactors = FALSE)
    rownames(met) <- met$drugID
    assign(esetNameP, new("ExpressionSet", exprs=targetsMatrix, phenoData=new("AnnotatedDataFrame", data=geneInfo),
			  featureData=new("AnnotatedDataFrame", data=met)))
    esetToGct1.3(eset=get(esetNameP), gctFile=paste0(esetNameP, ".gct"), path_to_write=path_to_write)
    
    return(paste0(esetNameP, ".gct"))
    
}    
    
#     cmpInd <- which(librarySigIDlist[, 1] == "LIB_5")
#     targetList <- NULL
# ##    load("/opt/raid10/genomics/data/lincs/proteinTargets/stitchLincs.RData")
# ##    This is from 2016
#     data(stitchLincs)
#     targetList <- split(stitchLincs[, "GeneID"], stitchLincs[, "lincsId"])
#     mtar <- match(alllincsPertID[cmpInd], names(targetList))
#     mtar1 <- mtar[!is.na(mtar)]
#     if (length(mtar1) > 0) {
#       if(!is.null(sigFile)) allProfileLabels <- allProfileLabels[-1]
#       allTargets <- unique(unlist(targetList[mtar1]))
#       # table where rows are profiles and columns are allTargets.
#       rows <- length(mtar1)
#       cols <- length(allTargets)
#       proteinTargetMatix <- matrix(c(rep(0, (rows * cols))), rows, cols)
#       for (i in 1:length(targetList[mtar1])) {
#         mIndex <- match(targetList[mtar1][[i]], allTargets)
#         proteinTargetMatix[i, mIndex] <- 10
#       }
#       colnames(proteinTargetMatix) <- paste(as.character(allTargets), unlist(AnnotationDbi::mget(as.character(allTargets), org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound = NA)), unlist(AnnotationDbi::mget(as.character(allTargets), org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound = NA)), sep = "::")
# #      colnames(proteinTargetMatix) <- paste(as.character(allTargets), unlist(mget(as.character(allTargets), getFromNamespace("org.Hs.egSYMBOL", ns="org.Hs.eg.db"), ifnotfound = NA)), unlist(mget(as.character(allTargets), org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound = NA)), sep = ":")
#       rownames(proteinTargetMatix) <- allProfileLabels[which(!is.na(mtar))]
#       # output this in treeview
#       if (rows == 1) {
#         rowClust <- fastcluster::hclust(dist(1:3))
#         nGenes <- rows
#         rowClust$merge <- rbind(c(-1, -2), matrix(c(seq(-3, -nGenes), seq(1:(nGenes - 2))), byrow = F, ncol = 2))
#         rowClust$height <- rep(1, nGenes - 1)
#         rowClust$order <- 1:nGenes
#         
#         colClust <- fastcluster::hclust(dist(1:3))
#         nSamples <- cols
#         colClust$merge <- rbind(c(-1, -2), matrix(c(seq(-3, -nSamples), seq(1:(nSamples - 2))), byrow = F, ncol = 2))
#         colClust$height <- rep(1, nSamples - 1)
#         colClust$order <- 1:nSamples
#         cluster = "none"
#       } else {
#         cluster <- "both"
#         if(ncol(proteinTargetMatix) == 1) {
# 	  cluster <- "row"
# 	  prot <- cbind(proteinTargetMatix, proteinTargetMatix)
# 	  colnames(prot)[2] <- " "
# # 	  prot<<-proteinTargetMatix
# 	  colClust <- fastcluster::hclust(dist(data.matrix(t(prot))), method = "complete")
# 	  rowClust <- fastcluster::hclust(dist(data.matrix((prot))), method = "complete")
# 	} else {
# 	  colClust <- fastcluster::hclust(dist(data.matrix(t(proteinTargetMatix))), method = "complete")
# 	  rowClust <- fastcluster::hclust(dist(data.matrix((proteinTargetMatix))), method = "complete")
# 	}
#       }
#       
#       # r2gtr(rowClust,file=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.gtr',sep=''),dec='.')
#       # r2atr(colClust,file=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.atr',sep=''),dec='.')
#       # r2cdt(hr=rowClust,hc=colClust,data=data.frame(rownames(proteinTargetMatix),rownames(proteinTargetMatix),proteinTargetMatix),labels=T,description=T,file=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.cdt',sep=''),dec = '.')
#       # call.treeview(paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),sep=''),ilincs=T)
#       # junk<-createJnlpFileGenomics(ilincs=T,fileName=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.jnlp',sep=''),location='/opt/raid10/genomics/Web/GenomicsPortals/ilincs', webLocation='http://eh3.uc.edu/genomics/GenomicsPortals/ilincs', targetFilesLocation='/opt/raid10/genomics/Web/GenomicsPortals/ilincs/treeviewFiles', targetFilesName=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.cdt',sep=''), targetWebLocation='http://eh3.uc.edu/genomics/GenomicsPortals/ilincs/treeviewFiles')
#       create_treeview_files(paste(sessionID, "_ProteinTargets", sep = ""), data.frame(rownames(proteinTargetMatix), rownames(proteinTargetMatix), proteinTargetMatix), 
# 			    ifCluster = cluster, isCentered = "FALSE", colClust, rowClust, path_to_write, sigDB=servSet$sigDB, chost=servSet$host, debuging=debuging)
#       esetNameP <- paste0(sessionID, "_ProteinTargets")
# #       met <- data.frame(meta=colnames(proteinTargetMatix), stringsAsFactors=FALSE)
#       met <- as.data.frame(t(apply(as.matrix(colnames(proteinTargetMatix)), 1, function(j) unlist(strsplit(j, "::")))), stringsAsFactors=F)
#       colnames(met) <- c("ID_geneid", "Name_GeneSymbol", "Description")
#       rownames(met) <- colnames(proteinTargetMatix) # <- met$geneSymbol # colnames(proteinTargetMatix)
#       
# #       fet <- data.frame(meta2 = rownames(proteinTargetMatix), stringsAsFactors=FALSE)
#       if(rownames(proteinTargetMatix)[1] == "uploadedSig") rownames(proteinTargetMatix)[1] <- paste(c("uploadedSig", NA,NA,NA,NA,NA),collapse="::")
#       fet <- as.data.frame(t(apply(as.matrix(rownames(proteinTargetMatix)), 1, function(j) unlist(strsplit(j, "::")))), stringsAsFactors=F)
#       fet <- fet[,-6]
#       colnames(fet) <- c("signatureID", "Perturbagen", "Concentration", "CellLine", "Time")
#       rownames(fet) <- rownames(proteinTargetMatix)
#       rownames(proteinTargetMatix) <- rownames(fet)
#       
#       assign(esetNameP, new("ExpressionSet", exprs=as.matrix(proteinTargetMatix), phenoData=new("AnnotatedDataFrame", data=met),
# 			    featureData=new("AnnotatedDataFrame", data=fet)))
#       ##save(list=esetNameP, file=paste0(path_to_write, esetNameP, ".RData"))
#       esetToGct1.3(eset=get(esetNameP), gctFile=paste0(esetNameP, ".gct"), path_to_write=path_to_write)
#       
#     } else {
#       # if protein target info present
#       # return error that none of the targets had protein target info
#       proteinRemark <- "not_generated_NoTargets"
#     }
# 
#     
#   ## return all the 3 results back to Java
#   res <- NULL
#   res[[1]] <- paste0(Remark, " Done")
#   names(res)[length(res)] <- "Remark"
#   res[[2]] <- sessionID
#   names(res)[length(res)] <- "SessionID"
#   if (is.null(proteinRemark)) {
#     res[[3]] <- paste("TreeView", sessionID, "_ProteinTargets", sep = "")
#   } else {
#     res[[3]] <- proteinRemark
#   }
#   names(res)[length(res)] <- "TreeViewProteinTargets"
#   if (is.null(onlyOneProfileRemark)) {
#     res[[4]] <- paste("TreeView", sessionID, "_ClusterTopGenes", sep = "")
#   } else {
#     res[[4]] <- onlyOneProfileRemark
#   }
#   names(res)[length(res)] <- "TreeViewClusterTopGenes"
#   res[[5]] <- paste("CorMatrix_", sessionID, ".xls", sep = "")
#   names(res)[length(res)] <- "GroupCorMatrix"
#   
#   if (is.null(topGeneExceedRemark)) {
#     res[[6]] <- ""
#   } else {
#     res[[6]] <- topGeneExceedRemark
#   }
#   names(res)[length(res)] <- "UnionExceed"
#   
#   if (is.null(NoOfGenesForHeatmapRemark)) {
#     res[[7]] <- ""
#   } else {
#     res[[7]] <- NoOfGenesForHeatmapRemark
#   }
#   names(res)[length(res)] <- "CappedOrInvalid"
#   
#   if (is.null(NoOfGenesForHeatmap)) {
#     res[[8]] <- ""
#   } else {
#     res[[8]] <- NoOfGenesForHeatmap
#   }
#   names(res)[length(res)] <- "UserEnteredGenes"
#   
#   res[[9]] <- length(topGenes)
#   names(res)[length(res)] <- "UnionLength"
#   
#   res[[10]] <- 50
#   names(res)[length(res)] <- "DefaultNoOfGenesThreshold"
#   
#   res[[11]] <- 1000
#   names(res)[length(res)] <- "MaximumNoOfGenesThreshold"
#   print(paste(" *** ", onlyOneProfileRemark, sep = ""))
#   print(paste(" *** ", NoOfGenesForHeatmapRemark, sep = ""))
#   print(paste(" *** ", topGeneExceedRemark, sep = ""))
#   
#   return(res)
# }
