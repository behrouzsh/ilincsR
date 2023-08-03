
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

master_GroupLincsAnalysis <- function(idList, libName="LIB_5", path_to_write, queryID=NULL, 
					NoOfGenesForHeatmap=NULL, glist=NULL, debuging=FALSE, sigDB="ilincs_sigs", sigFile=NULL) {
#libName = "LIB_5"
##  library(org.Hs.eg.db)
  # library(org.Mm.eg.db)
  # library(org.Rn.eg.db)
  print(idList)
#  print(libName)
  print(NoOfGenesForHeatmap)
  # generate sessionID
#   sessionID <- paste(c(strsplit(date(), split = " ")[[1]], as.integer(runif(1) * 1e+07)), sep = "", collapse = "_")
  ## Mehdi: changed the sessionID
  sID <- unlist(strsplit(date(), " "))
  sessionID <- paste(c(sID[1:3], sID[5], sID[4], as.integer(runif(1) * 1e+07)), sep = "", collapse = "_") 
  sessionID <- gsub(pattern = ":", replacement = "_", x = sessionID)
  # parse idList
  idListProfiles <- unlist(strsplit(idList, ","))
  # open connection
  servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
  mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = servSet$sigDB, port = servSet$port, host = servSet$host, password = "public")
  # get all the library types based on idList
  sql <- paste("select libraryID,signatureID from signaturesMeta3 where signatureID IN ('", gsub(",", "','", idList), "')", sep = "")
  mysql_result <- DBI::dbSendQuery(mycon, sql)
  librarySigIDlist <- DBI::fetch(mysql_result, n = -1)
  librarySigIDlist <- librarySigIDlist[match(idListProfiles, librarySigIDlist$signatureID), ]
  Remark <- "All signatures found - "
  if (dim(librarySigIDlist)[1] != dim(na.omit(librarySigIDlist))[1]) Remark <- "Some signatures not found - "
  
  librarySigIDlist <- na.omit(librarySigIDlist)
  libraryList <- unique(librarySigIDlist$libraryID)
  ## libraryList <- split(librarySigIDlist[,2],librarySigIDlist[,1])
  ## get profiles info from database based on libraryType
  ## sql<-paste("select pvaluesTablePath,diffexpTablePath,probsTablePath,scoresTablePath from signatureLibraries where libraryID = '",libName,"'",sep="")
  libraryList_paste <- paste(libraryList, collapse = ",")
  libraryList_sql <- gsub(",", "','", libraryList_paste)
  sql <- paste("select libraryID,pvaluesTablePath,diffexpTablePath,probsTablePath,scoresTablePath from signatureLibraries where libraryID IN ('", libraryList_sql, "')", sep = "")
  rs <- DBI::dbSendQuery(mycon, sql)
  allPaths <- DBI::fetch(rs, n = -1)
  allPaths <- allPaths[match(libraryList, allPaths$libraryID), ]
#@  allPaths <- allPaths
#@  libraryList <- libraryList
  # load all profiles
  refProfiles <- NULL
  refPvalProbProfiles <- NULL
  allProfileLabels <- NULL
  alllincsPertID <- NULL
  for (p in 1:length(libraryList)) {
    libraryPath <- allPaths[p, "pvaluesTablePath"]
    libraryPathDiffExpr <- allPaths[p, "diffexpTablePath"]
    libraryProb <- allPaths[p, "probsTablePath"]
    libraryScore <- allPaths[p, "scoresTablePath"]
    # also obtain signatureMetadata
    # sql<-paste("select * from signaturesMeta where signatureID IN ('",gsub(",","','",idList),"')",sep="")
    pInd <- which(librarySigIDlist$libraryID == libraryList[p])
    profileList <- librarySigIDlist$signatureID[pInd]
    profileList <- unlist(paste(profileList, collapse = ","))
    ## profileList <- unlist(paste(unlist(libraryList[p]),collapse=","))
    profileList_sql <- gsub(",", "','", profileList)
    sql <- paste("select * from signaturesMeta3 where signatureID IN ('", profileList_sql, "')", sep = "")
    rs <- DBI::dbSendQuery(mycon, sql)
    data_table <- DBI::fetch(rs, n = -1)
    sm <- match(librarySigIDlist$signatureID[pInd], data_table$signatureID)
    data_table <- data_table[sm, ]
    lincsPertID <- data_table$lincsPertID
    libName <- libraryList[p]
    
    
    profileLabels <- apply(data_table[, c("signatureID","datasetID","antibodyTarget","compound", "concentration", "cellLine", "time", "treatment", "factor", "Level1", "Level2")], 1, function(x) paste(x, collapse = "::"))
    
    allProfileLabels <- c(allProfileLabels, profileLabels)
    alllincsPertID <- c(alllincsPertID, lincsPertID)
    # load Rdata object for the respective data
    ## Mehdi: for speeding up
#@    info <- list(idListProfiles=idListProfiles, libraryList=libraryList, allPaths=allPaths, libName=libName, 
#		    librarySigIDlist=librarySigIDlist, pInd=pInd, 
#		    allProfileLabels=allProfileLabels, alllincsPertID=alllincsPertID, profileList=profileList)
		    
  if (length(idListProfiles) >= 100) {
    if (libName == "LIB_3") {
      temprefProfiles <- load(libraryScore)
      temprefPvalProbProfiles <- load(libraryProb)
    } else {
      temprefProfiles <- load(libraryPathDiffExpr)
#@      test <- temprefProfiles
      temprefPvalProbProfiles <- load(libraryPath)
    }
    } else {
##@!	if (sigDB == "ilincs_new") {
	    sigTable <- switch(libraryList[p], 
				LIB_1 = "GdsSigData", 
				LIB_2 = "CmapSigData", 
				LIB_3 = "EncodeSigData", 
				LIB_5 = "LincsCpSigData", 
				LIB_6 = "LincsKdSigData", 
				LIB_7 = "LincsShSigData", 
				LIB_8 = "LincsRnaSeqSigData", 
				LIB_9 = "LincsP100SigData", 
				LIB_10 = "ctrsSigData",
				LIB_11 = "LincsOeSigData",
				LIB_12 = "DrugMatrixSigData",
				LIB_13 = "ebiSigData",
				LIB_14 = "pharmGenSigData")
  if (sigDB == "ilincsTest2") {
    sigTable <- sub("KdSigData", "KdSigData_Nov12", sigTable)
    sigTable <- sub("CpSigData", "CpSigData_Nov12", sigTable)
  }
#     siglib <- switch(libraryList[p], 
# 			LIB_1 = "GDS", LIB_2 = "CMAP", LIB_3 = "ENC", LIB_5 = "LINCSCP", LIB_6 = "LINCSKD", LIB_7 = "LINCSSH", 
# 			LIB_8 = c("EDS-1014", "LDS-1237", "LDS-1239"), LIB_9 = "LINCSTP", LIB_10 = "CTRS")
    siglib <- switch(libraryList[p], 
			LIB_1 = "GDS", LIB_2 = "CMAP", LIB_3 = "ENC", LIB_5 = "LINCSCP", LIB_6 = "LINCSKD", LIB_7 = "LINCSSH", 
			LIB_8 = "LINCSRS", LIB_9 = "LINCSTP", LIB_10 = "CTRS", LIB_11 = "LINCSOE", LIB_12 = "DM", LIB_13 = "EBI", LIB_14 = "PG")

print(paste(sigTable, siglib))
#     allids <- paste0(idListProfiles[grep(siglib, idListProfiles)], collapse = "','")
    ids <- idListProfiles[grep(siglib, idListProfiles)] # ids <- unique(sigData$signatureID)
    ids <- NULL
    for (lname in siglib) {
	ids <- c(ids , idListProfiles[grep(lname, idListProfiles)])
    }
    allids <- paste0(ids, collapse = "','")
    sigData <- DBI::dbGetQuery(mycon, paste("select * from ",sigTable," where signatureID in ('", allids,"');",sep=""))
#    sigData <- sigData[sigData$geneID != 16777215, ] ## This was for problem in database definitions for geneID column!!!
#@allids <- allids
#@ids <- ids
    ## start casting
    allgenes <- unique(sigData$geneID)
    sigData <- split(sigData, sigData$signatureID)
#@sigData <- sigData

    temprefProfiles <- data.frame(row.names=allgenes)
    temprefPvalProbProfiles <- data.frame(row.names=allgenes)
 	for (si in sigData) { #2:length(ids)) {
#@	si <- si
	m <- match(rownames(temprefProfiles), si$geneID)
	temprefProfiles <- cbind(temprefProfiles, si[m,3])
	colnames(temprefProfiles)[dim(temprefProfiles)[2]] <- si[1,1]
	
	temprefPvalProbProfiles <- cbind(temprefPvalProbProfiles, si[m, 4])
	colnames(temprefPvalProbProfiles)[dim(temprefPvalProbProfiles)[2]] <- si[1,1]
	}
  }

    
#    if (length(idListProfiles) > 1) {
#	for (s in 2:length(idListProfiles)) {
#	temprefProfiles <- cbind(temprefProfiles, DBI::dbGetQuery(mycon, paste("select * from ",sigTable," where signatureID='", idListProfiles[s],"';",sep=""))[,3])
#	colnames(temprefProfiles)[s] <- idListProfiles[s]
#	temprefPvalProbProfiles <- cbind(temprefPvalProbProfiles, DBI::dbGetQuery(mycon, paste("select * from ",sigTable," where signatureID='", idListProfiles[s],"';",sep=""))[,4])
#	colnames(temprefPvalProbProfiles)[s] <- idListProfiles[s]
#	}
 #   }
#  }
     ### end
     
  if (length(idListProfiles) >= 100) {
    temprefProfiles <- get(temprefProfiles)
    temprefPvalProbProfiles <- get(temprefPvalProbProfiles)
  }
#   temprefPvalProbProfiles<-temprefPvalProbProfiles
#   temprefProfiles<-temprefProfiles
   # subset profiles based on idList
##    ind <- match(librarySigIDlist[pInd, 2], colnames(get(temprefProfiles)))
    ind <- match(librarySigIDlist[pInd, 2], colnames(temprefProfiles))
    if (length(ind) == 0) {
      # none match
      # return with possible error message
      return("Matching signatures error!")
    }
##    temprefProfiles <- get(temprefProfiles)[, ind, drop = F]
##    temprefPvalProbProfiles <- get(temprefPvalProbProfiles)[, ind, drop = F]
    temprefProfiles <- temprefProfiles[, ind, drop = FALSE]
    temprefPvalProbProfiles <- temprefPvalProbProfiles[, ind, drop = FALSE]
#     temprefProfiles <- temprefProfiles[, !is.na(ind), drop = FALSE]
#     temprefPvalProbProfiles <- temprefPvalProbProfiles[, !is.na(ind), drop = FALSE]

    if (p == 1) {
      refProfiles <- temprefProfiles
      refPvalProbProfiles <- temprefPvalProbProfiles
    } else {
      # merge the two profiles by genes
#       refgenes <- rownames(refProfiles)
#       querygenes <- rownames(temprefProfiles)
#       mm <- match(refgenes, querygenes)
      mm <- match(rownames(refProfiles), rownames(temprefProfiles))
#       genes <- refgenes[which(!is.na(mm))]
      # temprefData <- temprefProfiles[mm[-(which(is.na(mm)))],]
      # temprefProbsData <- temprefPvalProbProfiles[mm[-(which(is.na(mm)))],]
      
#       refProfiles<-refProfiles
#       refPvalProbProfiles<-refPvalProbProfiles
#       temprefProfiles<-temprefProfiles
#       temprefPvalProbProfiles<-temprefPvalProbProfiles
      
      temprefData <- temprefProfiles[mm[!is.na(mm)], , drop=FALSE]
      temprefProbsData <- temprefPvalProbProfiles[mm[!is.na(mm)], ,drop=FALSE]
      
      refProfiles <- refProfiles[rownames(temprefData), , drop=FALSE]
      refPvalProbProfiles <- refPvalProbProfiles[rownames(temprefProbsData), , drop=FALSE]
      
      refProfiles <- cbind(refProfiles, temprefData)
      refPvalProbProfiles <- cbind(refPvalProbProfiles, temprefProbsData)
    }
  }  # end of "for loop" over multiple libraries
#@refProfiles <- refProfiles
########################################################
  DBI::dbDisconnect(mycon)
  
#   allProfileLabels <- sapply(allProfileLabels, function(j) gsub(" ", "_", j))
#   allProfileLabels <- sapply(allProfileLabels, function(j) gsub("-", "_", j))
  allProfileLabs <- as.data.frame(t(as.data.frame(strsplit(allProfileLabels, "::"))), stringsAsFactors=F)
  colnames(allProfileLabs) <- c("signatureID","datasetID","antibodyTarget","compound", "concentration", "cellLine", "time", "treatment", "factor", "Level1", "Level2")
  allProfileLabs <- allProfileLabs[, apply(allProfileLabs, 2, function(j) length(j[j=="NA"])!=length(j))]
  allProfileLabels <- apply(allProfileLabs, 1, function(j) paste(j, collapse="::"))

  if(!is.null(sigFile)) {
#     uploaded <- read.table(file = sigFile, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
    fl <- unlist(strsplit(sigFile, "/"))
    uploaded <- meanOverGene(uploadedSigProcess(fl[length(fl)], paste(fl[1:(length(fl)-1)], collapse="/"), write=FALSE, debuging=debuging)$Sig)
    rownames(uploaded) <- uploaded$geneID
    
    if (dim(uploaded)[2] == 4) {
      colnames(uploaded)[3:4] <- c("coefficients", "Pvals")
    }
    if (dim(uploaded)[2] == 3) {
      Remark <- paste0(Remark," Missing p-Values in uploaded signature")
#       res <- list(Remark=Remark)
#       return(res)
      colnames(uploaded)[3] <- "coefficients"
      uploaded$Pvals <- pnorm(-abs(uploaded$coefficients))
    }
    
#@!     mchSig <- na.omit(match(rownames(refProfiles), as.character(uploaded[,1])))
    mchSig <- intersect(rownames(refProfiles), rownames(uploaded))
    uploaded <- uploaded[mchSig, ]
#     refProfiles <- refProfiles[mchSig, ]
#@!     refProfiles <- cbind(uploadedSig=uploaded$coefficients, refProfiles[uploaded$geneID, , drop=FALSE])
    refProfiles <- cbind(uploadedSig=uploaded$coefficients, refProfiles[mchSig, , drop=FALSE])
#@!    refPvalProbProfiles <- cbind(uploadedSig=uploaded$Pvals, refPvalProbProfiles[uploaded$geneID, ])
    refPvalProbProfiles <- cbind(uploadedSig=uploaded$Pvals, refPvalProbProfiles[mchSig, , drop=FALSE])
    allProfileLabels <- c("uploadedSig", allProfileLabels)
    allProfileLabs <- rbind(rep("uploadedSig", ncol(allProfileLabs)), allProfileLabs)
    }
    if (dim(refProfiles)[1] <= 1) {
      Remark <- "Filtered 0/1 genes"
      res <- list(Remark=Remark)
      return(res)
    }
    
  
  onlyOneProfileRemark <- NULL
  if (length(idListProfiles) == 1  & is.null(sigFile)) { # !("uploadedSig" %in% colnames(refProfiles))
    # do not need to generate pair-wise excel OR treeview of concordant signatures.
    #onlyOneProfileRemark <- "only_one_profile"
    groupCorMatrix = as.matrix(1)
    rownames(groupCorMatrix) <- allProfileLabels
  } else {
    # Result1: compute pair-wise CC for the profiles
    # if encode, remove it for now before clustering
    if ("LIB_3" %in% libraryList) {
      if (length(libraryList) == 1) {
        groupCorMatrix <- cor(refProfiles, refProfiles, use = "pairwise.complete.obs")
        rownames(groupCorMatrix) <- allProfileLabels
      } else {
        encInd <- which(librarySigIDlist[, 1] == "LIB_3")
        groupCorMatrix <- cor(refProfiles[, -encInd], refProfiles[, -encInd], use = "pairwise.complete.obs")
        rownames(groupCorMatrix) <- allProfileLabels[-encInd]
      }
    } else {
      groupCorMatrix <- cor(refProfiles, refProfiles, use = "pairwise.complete.obs")
      rownames(groupCorMatrix) <- allProfileLabels
    }
  }
  groupCorMatrix <- rbind(rownames(groupCorMatrix), groupCorMatrix)
  write.table(groupCorMatrix, file = paste(path_to_write, "CorMatrix_", sessionID, ".xls", sep = ""), sep = "\t", row.names = T, col.names = F, append = TRUE, quote=TRUE)  # save excel file with correlations
  # Result2: get top genes from each profile and take union of all genes
  dataClust <- refProfiles
  dataClustPval <- refPvalProbProfiles
  topGeneExceedRemark <- NULL
  NoOfGenesForHeatmapRemark <- NULL
  
  ## Mehdi: added glist to function
 if (!is.null(glist)) topGenes <- unique(na.omit(unlist(strsplit(glist, ",")))) else {
  if (is.null(NoOfGenesForHeatmap)) {
    # it is the first iteration so user input not obtained
    NoofGenes <- 50
  } else if (is.na(as.numeric(NoOfGenesForHeatmap))) {
    # not number
    NoOfGenesForHeatmapRemark <- "InvalidNumberOfGenesFormat"
    NoofGenes <- 50  # set to default minimum 
  } else if (as.integer(NoOfGenesForHeatmap) != NoOfGenesForHeatmap) {
    # number is not integer
    NoOfGenesForHeatmapRemark <- "InvalidNumberOfGenesFormat"
    NoofGenes <- 50  # set to default minimum 
  } else if (as.numeric(NoOfGenesForHeatmap) > 1000) {
    # it is valid integer
    NoofGenes <- 1000  # cap it to 1000
    NoOfGenesForHeatmapRemark <- "NoOfGenesCapped"
  } else {
    # user entered a valid number of genes
    NoofGenes <- NoOfGenesForHeatmap
  }
  # if(libName== 'LIB_4') { NoofGenes <- NoOfGenesForHeatmap }
  
  if(!is.null(sigFile)) {
    topGenes <- getTopGenes(dataClustPval, NoofGenes, c("LIB_X", librarySigIDlist[, 1]))  # libList with no encode 
  } else {
    topGenes <- getTopGenes(dataClustPval, NoofGenes, librarySigIDlist[, 1])  # libList with no encode
  }
}
#   dataClustPval<-dataClustPval
#   NoofGenes<-NoofGenes
#   librarySigIDlist<-librarySigIDlist
#   topGenes<-topGenes
	#tmp1 <- org.Hs.eg.db::org.Hs.egSYMBOL
	#mapped_genes <- AnnotationDbi::mappedkeys(tmp1)
	#tmpSym <- as.list(tmp1[mapped_genes])
	#save(tmpSym, file="./ilincsR/data/geneID2Symbol.rda")
	
	#tmp2 <- org.Hs.eg.db::org.Hs.egGENENAME
	#mapped_genes <- AnnotationDbi::mappedkeys(tmp2)
	#tmpNam <- as.list(tmp2[mapped_genes])
	#save(tmpNam, file="./ilincsR/data/geneID2Name.rda")
	
#	load("/opt/raid10/genomics/mehdi/ilincs/gitHub/ilincsR/data/geneID2Symbol.rda")
#	load("/opt/raid10/genomics/mehdi/ilincs/gitHub/ilincsR/data/geneID2Name.rda")

#  topGenes <- topGenes[topGenes %in% names(tmpSym)]

topGenes <- topGenes[topGenes %in% rownames(dataClust)]
  if (length(topGenes) == 0) {
    NoOfGenesForHeatmapRemark <- 0
    topGeneExceedRemark <- "noCommon_Genes"
    onlyOneProfileRemark <- "noCommon_Genes"
    
  } else if (length(topGenes) <= 10000 & length(topGenes) != 0) {
    dataClust <- dataClust[topGenes, , drop = F]
#    dataClust <- dataClust[!is.na(match(rownames(dataClust), topGenes)), ]

    if (length(libraryList) > 1) {
      # heterogenous list
      # if encode, remove it for now before clustering
      if ("LIB_3" %in% libraryList) {
        encInd <- which(librarySigIDlist[, 1] == "LIB_3")
        dataClust <- dataClust[, -encInd, drop = F]
        allProfileLabels <- allProfileLabels[-encInd]
      }
    }
    # annotatons
    colnames(dataClust) <- allProfileLabels
##    rownames(dataClust) <- paste(topGenes, as.character(unlist(tmpSym[topGenes])), as.character(unlist(tmpNam[topGenes])), sep = ":")
    if (is.null(glist)) { 
    rownames(dataClust) <- paste(topGenes, unlist(AnnotationDbi::mget(topGenes, org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound=NA)), unlist(AnnotationDbi::mget(topGenes, org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound=NA)), sep = "::")
    } else {
    rownames(dataClust) <- paste(rownames(dataClust), unlist(AnnotationDbi::mget(rownames(dataClust), org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound=NA)), unlist(AnnotationDbi::mget(rownames(dataClust), org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound=NA)), sep = "::")
#    rownames(dataClust) <- paste(rownames(dataClust), geneid2symbol(paste(rownames(dataClust), collapse=",")), unlist(AnnotationDbi::mget(rownames(dataClust), org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound=NA)), sep = ":")
   }
#    rownames(dataClust) <- paste(topGenes, unlist(mget(topGenes, getFromNamespace("org.Hs.egSYMBOL", ns="org.Hs.eg.db"), ifnotfound = NA)), unlist(mget(topGenes, getFromNamespace("org.Hs.egGENENAME", ns="org.Hs.eg.db"), ifnotfound = NA)), sep = ":")
    # dataClust[which(dataClust =='NaN')]<-NA
    co <- apply(dataClust, 2, function(x) sum(is.na(x)))
    ro <- apply(dataClust, 1, function(x) sum(is.na(x)))
    missingRows <- which(ro <= (dim(dataClust)[2])/2)
    missingCols <- which(co <= (dim(dataClust)[1])/2)
    missingSigs <- sapply(colnames(dataClust), function(j) unlist(strsplit(j, "::"))[1])[-missingCols]
    Remark <- paste0(Remark, " signatures:", paste(missingSigs, collapse = ","), " had more than 50% missing data -")
    dataClust <- dataClust[missingRows, missingCols, drop = F]
    dataClust <- as.data.frame(dataClust)
    if (dim(dataClust)[2] == 1) {
      dataClust$Zero = 0
    }
    if ((dim(dataClust)[1]) > 1) {
      
      colClust <- fastcluster::hclust(dist(data.matrix(t(dataClust))), method = "complete")
      rowClust <- fastcluster::hclust(dist(data.matrix((dataClust))), method = "complete")
      
      # r2gtr(rowClust,file=paste(path_to_write,'TreeView',paste(sessionID,'_ClusterTopGenes',sep=''),'.gtr',sep=''),dec='.')
      # r2atr(colClust,file=paste(path_to_write,'TreeView',paste(sessionID,'_ClusterTopGenes',sep=''),'.atr',sep=''),dec='.')
      # r2cdt(hr=rowClust,hc=colClust,data=data.frame(rownames(dataClust),rownames(dataClust),dataClust),labels=T,description=T,file=paste(path_to_write,'TreeView',paste(sessionID,'_ClusterTopGenes',sep=''),'.cdt',sep=''),dec = '.')
      # call.treeview(paste(path_to_write,'TreeView',paste(sessionID,'_ClusterTopGenes',sep=''),sep=''),ilincs=T)
      # junk<-createJnlpFileGenomics(ilincs=T,fileName=paste(path_to_write,'TreeView',paste(sessionID,'_ClusterTopGenes',sep=''),'.jnlp',sep=''),location='/opt/raid10/genomics/Web/GenomicsPortals/ilincs', webLocation='http://eh3.uc.edu/genomics/GenomicsPortals/ilincs', targetFilesLocation='/opt/raid10/genomics/Web/GenomicsPortals/ilincs/treeviewFiles',targetFilesName=paste(path_to_write,'TreeView',paste(sessionID,'_ClusterTopGenes',sep=''),'.cdt',sep=''), targetWebLocation='http://eh3.uc.edu/genomics/GenomicsPortals/ilincs/treeviewFiles')
      
      create_treeview_files(paste(sessionID, "_ClusterTopGenes", sep = ""), data.frame(topGenes[as.vector(missingRows)], rownames(dataClust), dataClust), 
			    ifCluster = "both", isCentered = "FALSE", colClust, rowClust, path_to_write, sigDB=servSet$sigDB, chost=servSet$host, debuging=debuging)
#       tmp<-dataClust
      ## Mehdi: saved eset for heatmap
      esetName <- paste0("filteredeset_", sessionID)
      dataClust2 <- dataClust
      colnames(dataClust2) <- sapply(colnames(dataClust2), function(j) unlist(strsplit(j, "::"))[1])
      allProfileLabs <- allProfileLabs[allProfileLabs$signatureID %in% colnames(dataClust2), ]
      
#       phen <- as.data.frame(t(as.data.frame(strsplit(colnames(dataClust2), "::"))), stringsAsFactors=F)
      if(colnames(dataClust2)[2]=="Zero") {
	  rownames(allProfileLabs) <- colnames(dataClust2)[1]
	  dataClust2 <- dataClust2[,-2, drop=FALSE]
      } else {
          rownames(allProfileLabs) <- colnames(dataClust2)
      }
      
##@!      rownames(allProfileLabs) <- colnames(dataClust2) <- allProfileLabs$signatureID
#       colnames(phen) <- c("signatureID","datasetID","antibodyTarget","compound", "cellLine", "treatment", "factor", "Level1", "Level2")
#       phen <- phen[,apply(phen,2,function(j) all(j!="NA"))]
      fdata <- as.data.frame(t(apply(as.matrix(rownames(dataClust2)), 1, function(j) unlist(strsplit(j, "::")))), stringsAsFactors=F)
      colnames(fdata) <- c("ID_geneid", "Name_GeneSymbol", "DESCRIPTION")
      rownames(fdata) <- rownames(dataClust2)
#       fdata<-fdata
#       dataClust2<-dataClust2
      assign(esetName, new("ExpressionSet", exprs=as.matrix(dataClust2), phenoData=new("AnnotatedDataFrame", data=allProfileLabs),
			    featureData=new("AnnotatedDataFrame", data=fdata)))
      save(list=esetName, file=paste0(path_to_write, esetName, ".RData"))
      esetToGct1.3(eset=get(esetName), gctFile=paste0(esetName, ".gct"), path_to_write=path_to_write)
    }
    
##### for only one gene
    if ((dim(dataClust)[1]) == 1) {
      dataClust <- rbind(dataClust, dataClust)
#       dat <- dataClust
      colClust <- fastcluster::hclust(dist(data.matrix(t(dataClust))), method = "complete")
      rowClust <- fastcluster::hclust(dist(data.matrix((dataClust))), method = "complete")
      
      create_treeview_files(paste(sessionID, "_ClusterTopGenes", sep = ""), data.frame(topGenes[as.vector(missingRows)], rownames(dataClust), dataClust), 
			    ifCluster = "both", isCentered = "FALSE", colClust, rowClust, path_to_write, sigDB=servSet$sigDB, chost=servSet$host, debuging=debuging)
      esetName <- paste0("filteredeset_", sessionID)
      dataClust2 <- dataClust
      
      rownames(allProfileLabs) <- colnames(dataClust2)
      fdata <- as.data.frame(t(apply(as.matrix(rownames(dataClust2)), 1, function(j) unlist(strsplit(j, "::")))), stringsAsFactors=F)
      colnames(fdata) <- c("ID_geneid", "Name_GeneSymbol", "DESCRIPTION")
      rownames(fdata) <- rownames(dataClust2)
      rownames(dataClust2)[2] <- rownames(fdata)[2] <- " "
      assign(esetName, new("ExpressionSet", exprs=as.matrix(dataClust2[1,,drop=F]), phenoData=new("AnnotatedDataFrame", data=allProfileLabs),
			    featureData=new("AnnotatedDataFrame", data=fdata[1,,drop=F])))
      save(list=esetName, file=paste0(path_to_write, esetName, ".RData"))
      esetToGct1.3(eset=get(esetName), gctFile=paste0(esetName, ".gct"), path_to_write=path_to_write)
    }
#####

  } else {
    # union of specified number of top genes from each profile exceeds certian limit
    # onlyOneProfileRemark <- paste(NoOfGenesForHeatmap,' number of top genes fetched ', length(topGenes),'genes', sep='')
    topGeneExceedRemark <- "topGenesExceeded"
    onlyOneProfileRemark <- "topGenesExceeded"
    
  }  # more than 1 profiles
# #   dataClust<-dataClust
# #   dataClustPval<-dataClustPval
  
  # Result 3: combine all targets (union)
  # if libType is ilincs and profile is compound, get targets info
  proteinRemark <- NULL
  if (length(grep("LIB_5", libraryList)) == 1) {
    # if(libName== 'LIB_5'){
    cmpInd <- which(librarySigIDlist[, 1] == "LIB_5")
    targetList <- NULL
##    load("/opt/raid10/genomics/data/lincs/proteinTargets/stitchLincs.RData")
##    This is from 2016
    data(stitchLincs)
    targetList <- split(stitchLincs[, "GeneID"], stitchLincs[, "lincsId"])
    mtar <- match(alllincsPertID[cmpInd], names(targetList))
    mtar1 <- mtar[!is.na(mtar)]
    if (length(mtar1) > 0) {
      if(!is.null(sigFile)) allProfileLabels <- allProfileLabels[-1]
      allTargets <- unique(unlist(targetList[mtar1]))
      # table where rows are profiles and columns are allTargets.
      rows <- length(mtar1)
      cols <- length(allTargets)
      proteinTargetMatix <- matrix(c(rep(0, (rows * cols))), rows, cols)
      for (i in 1:length(targetList[mtar1])) {
        mIndex <- match(targetList[mtar1][[i]], allTargets)
        proteinTargetMatix[i, mIndex] <- 10
      }
      colnames(proteinTargetMatix) <- paste(as.character(allTargets), unlist(AnnotationDbi::mget(as.character(allTargets), org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound = NA)), unlist(AnnotationDbi::mget(as.character(allTargets), org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound = NA)), sep = "::")
#      colnames(proteinTargetMatix) <- paste(as.character(allTargets), unlist(mget(as.character(allTargets), getFromNamespace("org.Hs.egSYMBOL", ns="org.Hs.eg.db"), ifnotfound = NA)), unlist(mget(as.character(allTargets), org.Hs.eg.db::org.Hs.egGENENAME, ifnotfound = NA)), sep = ":")
      rownames(proteinTargetMatix) <- allProfileLabels[which(!is.na(mtar))]
      # output this in treeview
      if (rows == 1) {
        rowClust <- fastcluster::hclust(dist(1:3))
        nGenes <- rows
        rowClust$merge <- rbind(c(-1, -2), matrix(c(seq(-3, -nGenes), seq(1:(nGenes - 2))), byrow = F, ncol = 2))
        rowClust$height <- rep(1, nGenes - 1)
        rowClust$order <- 1:nGenes
        
        colClust <- fastcluster::hclust(dist(1:3))
        nSamples <- cols
        colClust$merge <- rbind(c(-1, -2), matrix(c(seq(-3, -nSamples), seq(1:(nSamples - 2))), byrow = F, ncol = 2))
        colClust$height <- rep(1, nSamples - 1)
        colClust$order <- 1:nSamples
        cluster = "none"
      } else {
        cluster <- "both"
        if(ncol(proteinTargetMatix) == 1) {
	  cluster <- "row"
	  prot <- cbind(proteinTargetMatix, proteinTargetMatix)
	  colnames(prot)[2] <- " "
# 	  prot<<-proteinTargetMatix
	  colClust <- fastcluster::hclust(dist(data.matrix(t(prot))), method = "complete")
	  rowClust <- fastcluster::hclust(dist(data.matrix((prot))), method = "complete")
	} else {
	  colClust <- fastcluster::hclust(dist(data.matrix(t(proteinTargetMatix))), method = "complete")
	  rowClust <- fastcluster::hclust(dist(data.matrix((proteinTargetMatix))), method = "complete")
	}
      }
      
      # r2gtr(rowClust,file=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.gtr',sep=''),dec='.')
      # r2atr(colClust,file=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.atr',sep=''),dec='.')
      # r2cdt(hr=rowClust,hc=colClust,data=data.frame(rownames(proteinTargetMatix),rownames(proteinTargetMatix),proteinTargetMatix),labels=T,description=T,file=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.cdt',sep=''),dec = '.')
      # call.treeview(paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),sep=''),ilincs=T)
      # junk<-createJnlpFileGenomics(ilincs=T,fileName=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.jnlp',sep=''),location='/opt/raid10/genomics/Web/GenomicsPortals/ilincs', webLocation='http://eh3.uc.edu/genomics/GenomicsPortals/ilincs', targetFilesLocation='/opt/raid10/genomics/Web/GenomicsPortals/ilincs/treeviewFiles', targetFilesName=paste(path_to_write,'TreeView',paste(sessionID,'_ProteinTargets',sep=''),'.cdt',sep=''), targetWebLocation='http://eh3.uc.edu/genomics/GenomicsPortals/ilincs/treeviewFiles')
      create_treeview_files(paste(sessionID, "_ProteinTargets", sep = ""), data.frame(rownames(proteinTargetMatix), rownames(proteinTargetMatix), proteinTargetMatix), 
			    ifCluster = cluster, isCentered = "FALSE", colClust, rowClust, path_to_write, sigDB=servSet$sigDB, chost=servSet$host, debuging=debuging)
      esetNameP <- paste0(sessionID, "_ProteinTargets")
#       met <- data.frame(meta=colnames(proteinTargetMatix), stringsAsFactors=FALSE)
      met <- as.data.frame(t(apply(as.matrix(colnames(proteinTargetMatix)), 1, function(j) unlist(strsplit(j, "::")))), stringsAsFactors=F)
      colnames(met) <- c("ID_geneid", "Name_GeneSymbol", "Description")
      rownames(met) <- colnames(proteinTargetMatix) # <- met$geneSymbol # colnames(proteinTargetMatix)
      
#       fet <- data.frame(meta2 = rownames(proteinTargetMatix), stringsAsFactors=FALSE)
      if(rownames(proteinTargetMatix)[1] == "uploadedSig") rownames(proteinTargetMatix)[1] <- paste(c("uploadedSig", NA,NA,NA,NA,NA),collapse="::")
      fet <- as.data.frame(t(apply(as.matrix(rownames(proteinTargetMatix)), 1, function(j) unlist(strsplit(j, "::")))), stringsAsFactors=F)
      fet <- fet[,-6]
      colnames(fet) <- c("signatureID", "Perturbagen", "Concentration", "CellLine", "Time")
      rownames(fet) <- rownames(proteinTargetMatix)
      rownames(proteinTargetMatix) <- rownames(fet)
      
      assign(esetNameP, new("ExpressionSet", exprs=as.matrix(proteinTargetMatix), phenoData=new("AnnotatedDataFrame", data=met),
			    featureData=new("AnnotatedDataFrame", data=fet)))
      ##save(list=esetNameP, file=paste0(path_to_write, esetNameP, ".RData"))
      esetToGct1.3(eset=get(esetNameP), gctFile=paste0(esetNameP, ".gct"), path_to_write=path_to_write)
      
    } else {
      # if protein target info present
      # return error that none of the targets had protein target info
      proteinRemark <- "not_generated_NoTargets"
    }
  } else {
    proteinRemark <- "not_generated_WrongLibrary"
  }
  
  ## return all the 3 results back to Java
  res <- NULL
  res[[1]] <- paste0(Remark, " Done")
  names(res)[length(res)] <- "Remark"
  res[[2]] <- sessionID
  names(res)[length(res)] <- "SessionID"
  if (is.null(proteinRemark)) {
    res[[3]] <- paste("TreeView", sessionID, "_ProteinTargets", sep = "")
  } else {
    res[[3]] <- proteinRemark
  }
  names(res)[length(res)] <- "TreeViewProteinTargets"
  if (is.null(onlyOneProfileRemark)) {
    res[[4]] <- paste("TreeView", sessionID, "_ClusterTopGenes", sep = "")
  } else {
    res[[4]] <- onlyOneProfileRemark
  }
  names(res)[length(res)] <- "TreeViewClusterTopGenes"
  res[[5]] <- paste("CorMatrix_", sessionID, ".xls", sep = "")
  names(res)[length(res)] <- "GroupCorMatrix"
  
  if (is.null(topGeneExceedRemark)) {
    res[[6]] <- ""
  } else {
    res[[6]] <- topGeneExceedRemark
  }
  names(res)[length(res)] <- "UnionExceed"
  
  if (is.null(NoOfGenesForHeatmapRemark)) {
    res[[7]] <- ""
  } else {
    res[[7]] <- NoOfGenesForHeatmapRemark
  }
  names(res)[length(res)] <- "CappedOrInvalid"
  
  if (is.null(NoOfGenesForHeatmap)) {
    res[[8]] <- ""
  } else {
    res[[8]] <- NoOfGenesForHeatmap
  }
  names(res)[length(res)] <- "UserEnteredGenes"
  
  res[[9]] <- length(topGenes)
  names(res)[length(res)] <- "UnionLength"
  
  res[[10]] <- 50
  names(res)[length(res)] <- "DefaultNoOfGenesThreshold"
  
  res[[11]] <- 1000
  names(res)[length(res)] <- "MaximumNoOfGenesThreshold"
  print(paste(" *** ", onlyOneProfileRemark, sep = ""))
  print(paste(" *** ", NoOfGenesForHeatmapRemark, sep = ""))
  print(paste(" *** ", topGeneExceedRemark, sep = ""))
  
  return(res)
}
