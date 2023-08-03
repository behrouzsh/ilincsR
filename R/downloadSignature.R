
#' A function to download Signatures.
#'
#' This function allows you to download a list of Signatures from Genomics Portals database.
#' @param sigID This is a string object. A signature ID or a list of signature IDs separated by commas.
#' @param path_to_write The path to save the signature/s as a ".csv" file. The default is set to a temp folder: "/mnt/raid/tmp/"
#' @param display If user wants to get the signatures as an R object (data.frame) otherwise it will only be saved.
#' @param geneNames If user wants to get gene names along with Entrez IDs in the downloaded signature(s).
#' @param pValue Cut off based on p-value.
#' @param foldchange Cut off based on fold change.
#' @param sigDB Which signature version user wants to use in the analysis.
#' @param noOfGenes Number of top genes in the signature based on their p-values.
#' @param glist Subset the signature based on a gene list.
#' @param write If need to save the tab delimited file.
#' @param debuging This is for debuging purposes only. 
#' @param authors M. Fazel-Najafabadi
#' @keywords Signature
#' @export 
#' @examples
#'	## not run
#'	res <- downloadSignature(sigID="ENC_1,ENC_10,ENC_15", path_to_write="my/path/to/save/", display=TRUE, debuging=FALSE)
#'	## end not run

downloadSignature <- function(sigID, path_to_write="/mnt/raid/tmp/", display=TRUE, geneNames=TRUE, pValue=NULL, foldchange=NULL,
				sigDB="ilincs_sigs", species="Hs", noOfGenes="all", glist=NULL, write=TRUE, debuging=FALSE)
{
#  require(RMySQL)
res <- list()
if(!is.null(pValue)) pValue <- as.numeric(pValue)
if(!is.null(foldchange)) foldchange <- as.numeric(foldchange)

if (noOfGenes == 0) { res$fileName <- "No genes selected, No file generated"; res$signature <- NA; return(res) }
if (!is.null(glist)) glist <- unlist(strsplit(glist, ","))
fileName <- paste(c("sig", strsplit(date(),split=" ")[[1]], as.integer(runif(1)*10e6)),sep="",collapse="_")
# fileName <- paste0(gsub(":","_", fileName), ".xls")
fileName <- gsub(":","_", fileName)
res$fileName <- fileName

sigs <- sigID
sigID <- unlist(strsplit(sigID, ","))
if(length(grep("ENC_", sigID)) != 0 & length(grep("ENC_", sigID)) != length(sigID)) { res$fileName <- "Can not combine ENC with others!"; res$signature <- NA; return(res) }
#print(sigID)
signature <- NULL
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="10.165.4.231", port = 4040, password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')}

#   connSig <- DBI::dbConnect(RMySQL::MySQL(), host="192.168.1.102",user="public",password="public",dbname="ilincsTest2")

#@if (test) db <- "ilincs_new" else db <-"ilincsTest2"  ## needs to be fixed after switching new Sig libraries
servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
connSig <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = servSet$sigDB, port = servSet$port, host = servSet$host, password = "public")
#@connSig <- DBI::dbConnect(RMySQL::MySQL(), host="192.168.1.102",user="public",password="public",dbname=db)
#   connGene <- dbConnect(MySQL(), host="192.168.1.102",user="public",password="public",dbname="GeneDB")
for (sig in sigID) {
  sigType <- strsplit(sig,split="_")[[1]][1]
#  print(sigType)
  sigTable <- switch(sigType,
#		      LINCSKD = "LincsKdSigData_Nov12",
		      LINCSKD = "LincsKdSigData",
#		      LINCSCP = "LincsCpSigData_Nov12",
		      LINCSCP = "LincsCpSigData",
		      LINCSSH = "LincsShSigData",
		      `EDS-1014` = "LincsRnaSeqSigData",
		      `LDS-1237` = "LincsRnaSeqSigData",
		      `LDS-1239` = "LincsRnaSeqSigData",
		      LINCSRS = "LincsRnaSeqSigData",
		      GDS = "GdsSigData",
		      CMAP = "CmapSigData",
		      ENC = "EncodeSigData",
		      LINCSTP = "LincsP100SigData",
		      LINCSOE = "LincsOeSigData",
		      DM = "DrugMatrixSigData",
		      CTRS = "ctrsSigData",
		      EBI = "ebiSigData",
		      PG = "pharmGenSigData"
		      )
  
  cols <- switch(sigType,
		      LINCSKD = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      LINCSCP = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      LINCSSH = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      `EDS-1014` = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      `LDS-1237` = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      `LDS-1239` = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      LINCSRS = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      GDS = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      CMAP = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      ENC = "signatureID,geneID,score,prob,geneID",
		      LINCSTP = "signatureID,geneID,Diff_Exp,P_Value,probeID",
		      LINCSOE = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      DM = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      CTRS = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      EBI = "signatureID,geneID,Diff_Exp,P_Value,geneID",
		      PG = "signatureID,geneID,Diff_Exp,P_Value,geneID"
		      )

  if (sigDB=="ilincsTest2") {
    sigTable <- sub("KdSigData", "KdSigData_Nov12", sigTable)
    sigTable <- sub("CpSigData", "CpSigData_Nov12", sigTable)
  }
	tmp <- DBI::dbGetQuery(connSig, paste0("select ", cols, " from ",sigTable," where signatureID='",sig,"';"))
	colnames(tmp) <- c( "signatureID", "geneID", "Diff_Exp", "P_Value", "PROBE")
# 	if (sigType == "ENC") {
# 	  tmp <- tmp[, c("signatureID", "geneID", "score", "prob")]
# 	} else {
# 	  tmp <- tmp[, c("signatureID", "geneID", "Diff_Exp", "P_Value")]
# 	}
	    if (!is.null(glist)) {
			tmp <- tmp[tmp$geneID %in% glist, ]
	    }
	    if (!is.null(pValue)) {
			tmp <- tmp[!is.na(tmp[,4]), ]
			tmp <- tmp[tmp[,4] <= pValue, ]
	    }
	    if (!is.null(foldchange)) {
		if(length(foldchange)==1) {
		    tmp <- tmp[abs(tmp[,3]) >= foldchange, ]
		} else {
		    tmp <- tmp[tmp[,3] <= foldchange[1] | tmp[,3] >= foldchange[2], ]
		}
	    }
	    
	    
	    if (noOfGenes != "all") {
# 		    if (sigType == "ENC") {
# 			tmp <- tmp[order(tmp$prob),]
# 			tmp <- tmp[1:min(noOfGenes, dim(tmp)[1]),]
# 		    } else {
			tmp <- tmp[order(tmp$P_Value),]
			tmp <- tmp[1:min(noOfGenes, dim(tmp)[1]),]
# 		    }
	    }
	    if(sigType != "LINCSTP") tmp[,5] <- NA
	    
  signature <- rbind(signature, tmp)
#  print(sigTable)
}
## signature table colnames in database
# signatureID	geneID		score		prob
# signatureID	geneID		Diff_Exp	P_Value
# signatureID	geneID	probeID	Diff_Exp	P_Value

# if (noOfGenes == 0) { res$fileName <- "No genes selected, No file generated"; res$signature <- NA; return(res) }
# if (noOfGenes != "all") {
# 	signature <- signature[order(signature$signatureID, signature$P_Value),]#
# 	signature <- signature[1:min(noOfGenes, dim(signature)[1]),]
# 	
# }

#DBI::dbClearResult(RMySQL::dbListResults(connSig)[[1]])
DBI::dbDisconnect(connSig)
#   dbDisconnect(connGene)

if (dim(signature)[1]==0) { res$fileName <- "No genes found, No file generated"; res$signature <- NA; return(res) }
if (geneNames) {
    allids <- unique(signature$geneID)
    allnames <- geneid2symbol(paste(allids,collapse=","), species="Hs", debuging=debuging)
    mm <- match(signature$geneID, allnames$GeneID)
    signature$geneName <- allnames[mm, "Symbol"]
    signature <- signature[, c(1, 5, 2, 6, 3, 4)] 
# }
    if (sigType == "ENC") {
#  	colnames(signature) <- c("signatureID", "GeneID", "GeneNames", "Score", "Probability")
	colnames(signature) <- c("signatureID", "PROBE", "ID_geneid", "Name_GeneSymbol", "Score", "Probability")
	} else {
#  	colnames(signature) <- c("signatureID", "GeneID", "GeneNames", "coefficients", "Pvals")
	colnames(signature) <- c("signatureID", "PROBE", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")			
	}
} else {
    if (sigType == "ENC") {
#  	colnames(signature) <- c("signatureID", "GeneID", "GeneNames", "Score", "Probability")
	colnames(signature) <- c("signatureID", "PROBE", "ID_geneid", "Score", "Probability")
	} else {
#  	colnames(signature) <- c("signatureID", "GeneID", "GeneNames", "coefficients", "Pvals")
	colnames(signature) <- c("signatureID", "PROBE", "ID_geneid", "Value_LogDiffExp", "Significance_pvalue")			
	}
}
	
res$signature <- signature
if(write) {
      write.table(signature, file=paste0(path_to_write, res$fileName, ".xls"), quote=FALSE, sep="\t", row.names=F, col.names=T)
      
##      sigid <- sapply(sigID, function(j) unlist(strsplit(j, "_"))[1])
##      if (length(unique(sigid))==1) {
	  
	  if (sigType == "ENC") {
	  edat <- reshape2::recast(signature, ID_geneid~signatureID, measure.var='Score')
	  } else {
	  edat <- reshape2::recast(signature, ID_geneid~signatureID, measure.var='Value_LogDiffExp')
	  }

	  rownames(edat) <- edat$ID_geneid
	  edat <- edat[,-1, drop=FALSE]
	  pdat <-  getSignatureMeta(sigs, debuging=debuging, sigDB=sigDB)
	  pdat <- pdat[, !apply(pdat, 2, function(j) all(is.na(j)))]
	  rownames(pdat) <- pdat$signatureID
	  pdat <- pdat[, -1]
	  data(geneinfo)
	  mm <- match(rownames(edat), as.character(geneinfo$GeneID))
	  fdat <- data.frame(ID_geneid=rownames(edat), Name_GeneSymbol=geneinfo$Symbol[mm], description=geneinfo$description[mm], stringsAsFactors=FALSE)
	  rownames(fdat) <- fdat$ID_geneid
	  loadNamespace("Biobase")
	  assign("eset", new("ExpressionSet", exprs=as.matrix(edat), 
	      phenoData=new("AnnotatedDataFrame", data=pdat), 
	      featureData=new("AnnotatedDataFrame", data=fdat)))
# 	  esetToGct1.3(eset=eset, gctFile=sub(".xls", ".gct", res$fileName), path_to_write=path_to_write)
	  esetToGct1.3(eset=eset, gctFile=res$fileName, path_to_write=path_to_write)
##      }
    if (path_to_write == "/mnt/raid/tmp/") {
	print(paste0("Your signatures are saved in: [www/dev].ilincs.org/tmp/", res$fileName)) 
	} else {
	print(paste0("Your signatures are saved in: ",path_to_write, res$fileName))
	}
    if (display) return(res) else return(res$fileName)
} else {
    return(signature)
}

}
