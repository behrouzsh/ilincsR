
#' A function for visualizing, analysing and further investigating on a next-gene experiment.
#'
#' This function allows you to run multiple analysis and is one of the main work flows in GenomicsPortals.
#' @param exp This argument is the experiment name, defined as a string corresponding to GenomicsPortals experiment name. The default is set to NULL.
#' @param prop prop defines which property of the experiment is going to be used for statistical analysis to specify the groups of samples.
#' @param filterchk filterchk is in the form of a string. It can have multiple property of the pData to filter
#'	based on. Each pair of property:value should be separated by a comma ",". Pairs are saparated
#'	by colon ":". An example ot filterchk is :
#'	filterchk="property1:value1,property1:value2,property2:value1"
#' @param basement basement A two value vector containing "baseline" and "treatment".
#' @param includeORexclude This argument basically is designed to filter the pData based on what 
#'	filterchk is or selected property of eset when is set to "1", or the reverse selection of filterchk 
#'	when set to "2". 
#'	It also can have NULL or "n" values when there is no property selected, in this case 
#'	the function will return the original ExpressionSet.
#' @param ifCluster This parameter defines in which way the data should be clustered. The default value in this work flow is set to "Genes" which clusters only genes or rows in the data matrix.
#' @param up,down Specifies the maximum number of up and down regulated genes.
#' @param window.size ...
#' @param path_to_write This necessary parameter specifies the path which user wants to save the 
#'	output files (heatmaps, TreeView files and some other output files) to be saved in.
#' @param pORqVal Defines if "p" values are used in the analysis or "q" values.
#' @param libName Is a parameter to specify one of five Signature libraries in the GenomicsPortals database to calculate the concordance of generated signatures during the analysis against.
#' @param userUploadedProfilePath Is required if user wants to use his/her own provided gene Profiles in the analysis. It should be in the form of a string containig the path to Profile text file.
#'	This text file should contain geneIDs in the first column, Coefficients in the second and pValues in the third column (It's also possible to submit the file without pValues).
#' @param queryGenelist This is a string parameter consisting of user provided gene IDs separated by commas to be used in querying the dataset or experiment.
#' @param previousID This will be used if user wants to reanalyse a new dataset or experiment using the gene list derived from first round of analysis.
#' @param reAnalyse Should be set to TRUE only if same/another experiment or dataset is needed to be analysed based on the results of the first round of analysis.
#' @param write write argument is used (in downstreem function "parse_res") to save the SignatureTable for further uses.
#'	 It should be set to TRUE, as it is by default, to save the SignatureTable in the form of txt. 
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param authors M. Fazel-Najafabadi
#' @keywords exp, prop, reAnalyse.
#' @export 
#' @examples
#' ## not run
#' # For creating a signature:
#' res <- master_LincsDataAnalysis(exp="EDS-1013", prop="ER", filterchk="ER:+,ER:-", 
#'	basement=NULL, includeORexclude="1", ifCluster="both", pORqVal="p", 
#'	path_to_write="/patho/to/your/output/directory/", 
#'	userUploadedProfilePath=NULL, up=4000,	down=1000, window.size=50)
#'
#' # For creating a signature and calculating concordances to a library of signatures:
#' res <- master_LincsDataAnalysis(exp="EDS-1013", prop="ER", filterchk="ER:+,ER:-", 
#'	basement=NULL, includeORexclude="1", ifCluster="both", pORqVal="p", 
#'	path_to_write="/patho/to/your/output/directory/", libName="LIB_5",
#'	userUploadedProfilePath=NULL, up=4000,	down=1000, window.size=50)
#'
#' # For reanalysing:
#' res <- master_LincsDataAnalysis(exp="EDS-1013",prop="ER", filterchk= "ER:+,ER:-", 
#'	basement= NULL, includeORexclude="1", ifCluster="both", pORqVal="p", 
#'	path_to_write="/patho/to/your/output/directory/", libName="LIB_5",
#'	userUploadedProfilePath=NULL, up=4000,	down=1000, window.size=50, 
#'	previousID="find previousID from your saved results in the first step", reAnalyse=TRUE)
#'
#'	# previousID may look something like this: "Mon_Oct_26_23_49_56_2015_1084923"
#'	# "Date + a randomely generated number separated by under scores"
#' # For Random Set (RS) analysis pipeline:
#' glist <- "57535,11010,283987,2920,2113,10612,9619,3263,2300,1282,148327,84002,8038,114793,4059,6591,3576,51523,25803,8828,4478,
#'	26577,8537,8601,29940,80380,858,4067,718,10809,221061,7043,388,57222,8710,6288,3169,7494,123099,84959,5268,56938,301,
#'	5806,11094,8492,4015,10644,57057,7033,140876,5270,3624,8061,9734,7431,5292,7128,169792,10333,2921,8876,3908,399,10916,
#'	776,3767,5329,83451,1191,25837,199221,2902,6549,309,8991,222,55243,124220,10643,414318,23171,100131187,168544,22996,57134,
#'	2065,55612,26996,6566,2152,27092,54847,857,119,284119,1284,79917,3628,5087"
#' resRS <- master_LincsDataAnalysis(exp=NULL, libName="LIB_2", path_to_write= "/patho/to/your/output/directory/", queryGenelist=glist)
#'
#' # For uploading a signature and finding concordances to a library:
#' res <- master_LincsDataAnalysis(exp="sample.tsv",prop=NULL, filterchk= NULL, 
#'	basement= NULL, includeORexclude=NULL, ifCluster="both", pORqVal="p", 
#'	path_to_write="/patho/to/your/output/directory/", libName="LIB_5",
#'	userUploadedProfilePath="/patho/to/your/output/directory/", up=4000, down=1000, window.size=50)
#'
#' ## not run
#' @export 

master_LincsDataAnalysis <- function(exp=NULL,		prop=NULL,	filterchk="n",		basement=NULL,	includeORexclude=2,
				   ifCluster="Genes",	up=4000,	down=1000,		window.size=50,	path_to_write,
				   pORqVal="p",		libName=NULL,	userUploadedProfilePath=NULL, previousID=NULL, reAnalyse=FALSE, 
				   queryGenelist=NULL, pValue=NULL, foldchange=NULL, write=TRUE, 
				   homoloGenes=TRUE, debuging=FALSE, multigroup=FALSE, sigDB="ilincs_sigs", chost="ilincs.org", org="Hs") ## Why filters instead of filterStr? in javaScripts all filters
				   ## are called filterStr!
				   ## I changed ifCluster from none to both.
				   ## I set pORqVal to "p" for now!
				   # Mehdi: we need to changed the libName from "LIB_1" to NULL as defult too or add NULL as an option!
{

filterchk <- gsub("<NA>", "NA", filterchk)
MLinfo <- list(exp=exp, prop=prop, filterchk=filterchk, basement=basement, includeORexclude=includeORexclude, ifCluster=ifCluster,
		up=up, down=down, window.size=window.size, path_to_write=path_to_write, pORqVal=pORqVal, libName=libName, 
		userUploadedProfilePath=userUploadedProfilePath, queryGenelist=queryGenelist, previousID=previousID, reAnalyse=reAnalyse)
## Mehdi: log

#        print(paste("exp",exp))
#        print(paste("prop",prop))
#        print(paste("filterchk",filterchk))
#        print(paste("basement",basement))
#        print(paste("includeORexclude",includeORexclude))
#        print(paste("ifCluster",ifCluster))
#        print(paste("up",up))
#        print(paste("down",down))
#        print(paste("window.size",window.size))
#        print(paste("path_to_write",path_to_write))
#        print(paste("pORqVal",pORqVal))
#        print(paste("libName",libName))
#        print(paste("userUploadedProfilePath",userUploadedProfilePath))
#        print(paste("queryGenelist",queryGenelist))
#        print(paste("previousID",previousID))
#        print(paste("reAnalyse",reAnalyse))
#	 print(paste("before get pdata"))

##	library(CLEAN)
##	library(RMySQL)
##	library(gtools)
##	library(multicore)
##	library(parallel) ##@Mehdi changed the library multicore
##	library(Biobase)
  # mod_Rserve_prod.R is already sourced in diffGenesLincs_prod.R!!! - Nick
  # This "source" call will be instantly overwritten by the source in diffGenesLincs_prod! - Nick
  ## Mehdi: but not always! 
  #source("http://www.eh3.uc.edu/r/mod_Rserve_prod.R")            #why is is here instead web.xml
  #source("http://www.eh3.uc.edu/r/diffGenesLincs_prod.R")         #rtest for development, r for production
#  source("/opt/raid10/genomics/software/RPrograms/source/ilincs/dev/mod_Rserve_prod.R") # for testing - Mehdi
#  source("/opt/raid10/genomics/software/RPrograms/source/ilincs/dev/diffGenesLincs_prod.R") # for testing - Mehdi
##  source("/opt/raid10/genomics/mehdi/ilincs/code/R/iLincs_R_Scripts/mod_Rserve_prod.R") # for testing - Mehdi
##  source("/opt/raid10/genomics/mehdi/ilincs/code/R/iLincs_R_Scripts/diffGenesLincs_prod.R") # for testing - Mehdi

	# if exp=NULL, queryGenelist is should not be null. call genelist option.
        # else  if userUploadedProfilePath == NULL,  thne proceed with existing pipeline calling diffGenes.R function to prepare query profile
	#	else custom file upload option is active and "exp" will have name of the file and userUploadedProfilePath will have path of the file 
  ## Mehdi: recalling previous arguments
# # # # #   if(reAnalyse) {
# # # # #     load(paste(path_to_write, "MLinfo_filteredeset_", previousID, ".RData", sep=""))
# # # # #     exp <- MLinfo$exp
# # # # #     prop <- MLinfo$prop
# # # # #     filterchk <- MLinfo$filterchk
# # # # #     basement <- MLinfo$basement
# # # # #     includeORexclude <- MLinfo$includeORexclude
# # # # #     ifCluster <- MLinfo$ifCluster
# # # # #     up <- MLinfo$up
# # # # #     down <- MLinfo$down
# # # # #     window.size <- MLinfo$window.size
# # # # #     path_to_write <- MLinfo$path_to_write
# # # # #     pORqVal <- MLinfo$pORqVal
# # # # #     libName <- MLinfo$libName
# # # # #     userUploadedProfilePath <- MLinfo$userUploadedProfilePath
# # # # # #     queryGenelist <- MLinfo$queryGenelist
# # # # # #    previousID <- MLinfo$previousID
# # # # # #    reAnalyse <- MLinfo$reAnalyse
# # # # #   }
	
	if (!is.null(libName)) { ##Mehdi added
	
	print(paste("libName is ",libName,sep=""))
	if(is.null(exp)){
        	metacols <- switch(libName, 
			    LIB_1="datasetID,factor,Level1,Level2",
			    LIB_2="compound,concentration,cellLine",
			    LIB_3="treatment,antibodyTarget,cellLine",
			    LIB_5="compound,concentration,cellLine,time",
			    LIB_6="treatment,perturbagenID,cellLine,time",
			    LIB_7="treatment,perturbagenID,cellLine,time",
			    LIB_8="datasetID,factor,Level1,Level2",
			    LIB_10="perturbagenID,compound,cellLine",
			    LIB_11="treatment,perturbagenID,cellLine,time",
			    BING="compound,lincsPertID,concentration,cellLine,time"
			    )
		
		## call pipeline to conduct analysis for genelist
		queryGenelist <- unlist(strsplit(queryGenelist,","))
		usedGenes <- length(unique(queryGenelist))
		enrichTable <- computeGenelistEnrichment(queryGenelist, libName, debuging)
# 		enrichTable <- enrichTable
		if(dim(na.omit(enrichTable))[1] > 2){ remark <- "Done" 
		} else{ remark <- "Error in Random Set or No enriched genes found" 
		    enrichTable <- na.omit(enrichTable)}
		# Note: in other cases, sessionID is written by plot_with_properties_eset function in mod_Rserve_prod.R - Nick
		sessionID <- generateSID()
		
#@		tableNames <-paste(as.vector(enrichTable[,1]),collapse="DELIM") 
#@		zScores <-paste(as.vector(enrichTable[,3]),collapse="DELIM") 
#@                res <- data.frame(sessionID,remark,tableNames,zScores)
        	enrichTablePath <- paste(path_to_write,"enrichTable_",sessionID,".xls",sep="")
                sigScores <- data.frame(signature=enrichTable[, 1], score=round(enrichTable[, 3],4), stringsAsFactors=FALSE)
                
#@		colnames(res) <- c("sessionID","Remark","tableNames","zScores")
        	# save it as excel file and send the path
#@        	write.table(enrichTable[,-2],file=paste(path_to_write,"enrichTable_",as.character(res[[1]]),".xls",sep=""),sep="\t",row.names=F,col.names=F,append=T)
        	write.table(enrichTable[,-2],file=paste(path_to_write,"enrichTable_",sessionID,".xls",sep=""), sep="\t", quote=TRUE, row.names=F, col.names=F, append=T)
#@        	res[[(length(res)+1)]] <- paste(path_to_write,"enrichTable_",as.character(res[[1]]),".xls",sep="")
#@        	colnames(res)[dim(res)[2]] <- "corTablePath"
        	
	    sigmeta <- getSignatureMeta(prop=metacols, signatures=paste(sigScores$signature, collapse=","), debuging=debuging, sigDB=sigDB)
	    sigmatch <- match(sigScores$signature, sigmeta$signatureID)
	    sigmeta <- sigmeta[sigmatch,]
	    sigmeta <- sigmeta[,-1]
	    if (libName=="LIB_1" | libName=="LIB_8") {
	    sigScores <- cbind(sigScores, sigmeta[,1:4])
	    colnames(sigScores) <- c("tableNames", "zScores", "datasetID","factor","Level1","Level2")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
	    #colnames(sigScores)[(dim(sigScores)[2]-3):dim(sigScores)[2]] <- c("datasetID","factor","Level1","Level2")
	    } else if  (libName=="LIB_2") {
	    sigScores <- cbind(sigScores, sigmeta[,1:3])
	    colnames(sigScores) <- c("tableNames", "zScores", "compound","concentration","cellLine")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
	    #colnames(sigScores)[(dim(sigScores)[2]-2):dim(sigScores)[2]] <- c("compound","concentration","cellLine")
	    } else if (libName=="LIB_3") {
	    sigScores <- cbind(sigScores, sigmeta[,1:3])
	    colnames(sigScores) <- c("tableNames", "zScores", "treatment","antibodyTarget","cellLine")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
	    #colnames(sigScores)[(dim(sigScores)[2]-2):dim(sigScores)[2]] <- c("treatment","antibodyTarget","cellLine")
	    } else if (libName=="LIB_5" | libName=="BING") {
	    sigScores <- cbind(sigScores, sigmeta[,1:4])
	    colnames(sigScores) <- c("tableNames", "zScores", "compound","concentration","cellLine","time")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
	    #colnames(sigScores)[(dim(sigScores)[2]-3):dim(sigScores)[2]] <- c("compound","concentration","cellLine","time")
	    } else if (libName=="LIB_10") {
	    sigScores <- cbind(sigScores, sigmeta[,1:3])
	    colnames(sigScores) <- c("tableNames", "zScores", "perturbagenID", "compound", "cellLine")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
	    #colnames(sigScores)[(dim(sigScores)[2]-3):dim(sigScores)[2]] <- c("compound","concentration","cellLine","time")
	    } else {#	if (libName="LIB_6" | libName="LIB_7") {
	    sigScores <- cbind(sigScores, sigmeta[,1:4])
	    colnames(sigScores) <- c("tableNames", "zScores", "treatment","perturbagenID","cellLine","time")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
	    #sigScores[[(length(sigScores)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
	    #colnames(sigScores)[(dim(sigScores)[2]-3):dim(sigScores)[2]] <- c("treatment","perturbagenID","cellLine","time")
	    }
	res <- list(usedGenes=usedGenes , sessionID=sessionID, Remark=remark, enrichTablePath=enrichTablePath, sigScores=sigScores)
       	return(res)
	} else {
		if(!is.null(userUploadedProfilePath)){
# 			#upload custom profile from a file
# 			queryProfile1 <- try(read.table(file = paste(userUploadedProfilePath,"/", exp, sep=""), header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE))
# 			queryProfile2 <- try(read.table(file = paste(userUploadedProfilePath,"/", exp, sep=""), header=TRUE, quote="", sep=" ", stringsAsFactors=FALSE))
# 			queryProfile3 <- try(read.table(file = paste(userUploadedProfilePath,"/", exp, sep=""), header=TRUE, quote="", sep=";", stringsAsFactors=FALSE))
# 			if (dim(queryProfile2)[2] >= dim(queryProfile1)[2]) queryProfile <- queryProfile2 else queryProfile <- queryProfile1
# 			if (dim(queryProfile3)[2] >= dim(queryProfile)[2]) queryProfile <- queryProfile3
# 			if(dim(queryProfile)[1] == 0){
# 				res <- getDefaultErrorDataframe("empty file")
# 				return(res)
# 			}
# # 			if(colnames(queryProfile)[1:4] == c("ID", "PROBE", "SYMBOL", "DESCRIPTION")) queryProfile <- queryProfile[, c("ID", "Diff_Exp", "Pvalues")]
# 			if(identical(colnames(queryProfile)[1:4], c("ID", "PROBE", "SYMBOL", "DESCRIPTION"))) queryProfile <- queryProfile[, c("ID", "Diff_Exp", "Pvalues")]
# 			if(colnames(queryProfile)[1] == "signatureID") queryProfile <- queryProfile[,-1]
# 			if(dim(queryProfile)[2] == 4) colnames(queryProfile) <- c("geneID","geneName","coefficients","Pvals")
# 			if(dim(queryProfile)[2] == 4) queryProfile <- queryProfile[, c(1,3,4)]
# 			if(dim(queryProfile)[2] == 3) colnames(queryProfile) <- c("geneID","coefficients","Pvals")
# 			if(dim(queryProfile)[2] == 2) colnames(queryProfile) <- c("geneID","coefficients")
# 			if((dim(queryProfile)[2] %in% c(2,3)) & is.na(as.integer(queryProfile[1,1]))) {  ## sum(sapply(as.integer(queryProfile[,1]), is.na)) != dim(queryProfile)[1]
# 			    symid <- symbol2geneid(paste(queryProfile$geneID, collapse=","))
# 			    symid <- symid[!duplicated(symid$Symbol), , drop=FALSE]
# 			    si <- match(queryProfile$geneID, symid$Symbol)
# 			    queryProfile$geneID <- symid$GeneID[si]
# 			    queryProfile <- queryProfile[!is.na(queryProfile$geneID), ]
# 			    }
# 			if(dim(queryProfile)[2] != 2 & dim(queryProfile)[2] != 3 & dim(queryProfile)[2] != 4){
# 				res <- getDefaultErrorDataframe("wrong format")
# 				return(res)
# 			}
# 			queryProfile <- na.omit(queryProfile)
# 			#####Generating sessionID
# 			# Note: in other cases, sessionID is written by plot_with_properties_eset function in mod_Rserve_prod.R - Nick
# 			sessionID <- paste(c(strsplit(date(),split=" ")[[1]],as.integer(runif(1)*10e6)),sep="",collapse="_")
# 			sessionID <- gsub(":","_", sessionID)
# 			sessionID <- gsub("__","_", sessionID)
# 
# 			res <- getDefaultErrorDataframe("Done")  #dummy object
# 			res[1] <- sessionID
# 			res <- res[-(c(18,19,20))]
			if (length(grep("processedSig_", exp)) ==1) {
			  queryProfile0 <- try(read.table(file = paste(userUploadedProfilePath,"/", exp, sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)) # , quote=""
			} else {
			  queryProfile0 <- uploadedSigProcess(exp, userUploadedProfilePath, path_to_write=path_to_write, write=FALSE, debuging=debuging)$Sig
			}
			if (dim(queryProfile0)[2] == 3) {
			    colnames(queryProfile0) <- c("geneID","geneName","coefficients")
			} else {
			    colnames(queryProfile0) <- c("geneID","geneName","coefficients","Pvals")
			}
			queryProfile <- queryProfile0[, -2]
			sessionID <- generateSID()
			res <- getDefaultErrorDataframe("Done")  #dummy object
			res[1] <- sessionID
			res <- res[-(c(18,19,20))]
		} else{
	        	##@Mehdi Generating diff Gene list/profile
	        	res <- master_diffGenes(exp=exp, prop=prop, filterchk=filterchk, basement=basement, includeORexclude=includeORexclude,
						ifCluster=ifCluster, up=4000, down=1000, window.size=50, path_to_write=path_to_write,
						pORqVal=pORqVal, previousID=previousID, reAnalyse=reAnalyse, 
						glist=queryGenelist, pValue=pValue, foldchange=foldchange, MLinfo=MLinfo, homoloGenes=homoloGenes, 
						multigroup=multigroup, sigDB=sigDB, chost=chost, debuging=debuging)
	        	
			level1str <- res[[15]]
			##print(str(res)) # Mehdi 1 line
			if(as.character(res[[2]])== "NA"){ #error in anova or no genes found etc.
				res <- getDefaultErrorDataframe(res[[6]])  # remark from function
				res[[15]] <- level1str    #different if "no two levels"
				print(res)
				return(res)
			}
			queryProfile <- data.frame(geneID=unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[13])))," ")),
					coefficients=as.numeric(unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[12])))," "))), 
					Pvals=as.numeric(unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[17])))," "))), stringsAsFactors=FALSE)
		}
		#Merge multiple rows for the same geneIds
		if(dim(queryProfile)[2] == 3) unqIdPvals <- split(queryProfile[,"Pvals"], queryProfile[,"geneID"])
		if(dim(queryProfile)[2] == 3) mergePvals <- unlist(lapply(unqIdPvals, function(x) exp(mean(log(x), na.rm = TRUE)))) ## Mehdi: I changed it so we dont loose genes with p-val=NA
 		unqIdCoeffs <- split(queryProfile[,"coefficients"], queryProfile[,"geneID"])
		mergeCoeffs <- unlist(lapply(unqIdCoeffs, function(x) mean(as.numeric(x), na.rm = TRUE)))
		if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(geneID=names(unqIdCoeffs),coefficients=mergeCoeffs,Pvals=mergePvals)
                if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(geneID=names(unqIdCoeffs),coefficients=mergeCoeffs)
	}
	##@Mehdi Looking for the Human Homologus genes to use in analysis (if sp is not Hs)
	# get homologous Hs genes for query profile if specie is not Hs
# 	if(org != "Hs") {
	  homoloArr <- getHomologousGenes(as.character(queryProfile[,1]), org=org)
# 	  if(length(which(!is.na(homoloArr))) == 0){
# 		  res <- getDefaultErrorDataframe("No homolo genes")
# 	  }
	  if(length(!is.na(homoloArr[,3])) == 0){
		  res <- getDefaultErrorDataframe("No homolo genes")
	  }
	  queryProfile[,1] <- homoloArr[match(queryProfile[,1], homoloArr[,1]),3]
# 	}
	queryProfile <- queryProfile[which(!is.na(queryProfile[,1])),]
	rownames(queryProfile) <- queryProfile$geneID ## Mehdi: I added this later on

	# query database tables for the selected dblist to obtain the location of RData objects
	servSet <- getServerSettings(debuging=debuging, sigDB=sigDB)
#        mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=as.character(servSet[[1]]), host=as.character(servSet[[2]]), port = as.numeric(as.character(servSet[[3]])), password="public")
        mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=servSet$sigDB, host=servSet$host, port = servSet$port, password="public")
        sql <- paste("select diffexpTablePath,scoresTablePath,probsTablePath from signatureLibraries where libraryID = '",libName,"'",sep="")
        rs <- DBI::dbSendQuery(mycon, sql)
#        if (unlist(strsplit(rs, " "))[1] == "Error") {
#	    servSet <- getServerSettings()
#	    mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=as.character(servSet[[1]]), host=as.character(servSet[[2]]), port = as.numeric(as.character(servSet[[3]])), password="public")
#	    rs <- DBI::dbSendQuery(mycon, sql)
#        }
        dt <- DBI::fetch(rs, n = -1)
        libraryPath <- dt[1,"diffexpTablePath"]
	libraryScore <- dt[1,"scoresTablePath"]
	libraryProb <- dt[1,"probsTablePath"]
        DBI::dbDisconnect(mycon)

#	libraryPath <- "/opt/raid10/rdataNew/lincs/datasets/cmapDiffExp.RData"      ## obtain the path with Pvalue profiles 
	# load Rdata object for the respective data
	##@Mehdi loading the appropriate LIB
# # # # # 	if (libName=="LIB_3"){ 
# # # # # 	    refProfiles <- get(load(libraryScore)); refProbs <- get(load(libraryProb))
# # # # # 	} else { 
# # # # # 	    refProfiles <- get(load(libraryPath))
# # # # # 	}

	## New, getting data from ilincsRdata
	dat <- switch(libName, 
			    LIB_1="gdsDiffExp",
			    LIB_2="cmapDiffExp",
			    LIB_3=c("encodeScores", "encodeProbs"),
			    LIB_5="lincscpDiffExp",
			    LIB_6="lincskdDiffExp",
			    LIB_7="lincsshDiffExp",
			    LIB_8="lincsRNAseqDiffExp",
			    LIB_9="lincsP100DiffExp",
			    LIB_10="ctrsDiffExp"
		)
	      
	if (libName=="LIB_3"){ 
	    refProfiles <- get(dat[1], envir=.GlobalEnv); refProbs <- get(dat[2], envir=.GlobalEnv)
	} else { 
	    refProfiles <- get(dat, envir=.GlobalEnv)
	}
	####

	
	print(libraryPath)
	# match common genes for queyProfile and refProfiles
	refgenes <- rownames(refProfiles)
	querygenes <- as.character(queryProfile[,1]) ## Mehdi: why it was integer?!!!
# 	mm <- match(refgenes,querygenes)
# 	genes <- refgenes[which(!is.na(mm))]
# 	queryData <- queryProfile[mm[-(which(is.na(mm)))],]
	genes <- intersect(refgenes, querygenes)
	queryData <- queryProfile[genes, ]

	if (libName=="LIB_3" & dim(queryProfile)[2] == 3){
# 		refData <- refProfiles[which(!is.na(mm)), ,drop=FALSE] 
		refData <- refProfiles[genes, ,drop=FALSE] 
		no <- apply(refData,2,function(i) all(is.na(i)))
		refData <- refData[, !no ,drop=FALSE]
# 		refProbsData <- refProbs[which(!is.na(mm)), ,drop=FALSE]
		refProbsData <- refProbs[genes, ,drop=FALSE]
		refProbsData <- refProbsData[, !no ,drop=FALSE]
	} else {
# 		refData <- refProfiles[which(!is.na(mm)), ,drop=FALSE]
		refData <- refProfiles[genes, ,drop=FALSE]
		no <- apply(refData,2,function(i) all(is.na(i)))
		refData <- refData[, !no ,drop=FALSE]
	}

	if (dim(queryData)[1] == 0) {
		res <- getDefaultErrorDataframe("Filtered zero genes") # no overlapping genes between query and reference
		print("********************* Filtered zero genes!")
		return(res)
	}
	if (dim(queryData)[1] == 1) { ## Mehdi: it was only == 0 before!
		res$Remark <- paste0("Filtered one. Only geneID ", rownames(queryData), " is in common.") # no overlapping genes between query and reference
		print(paste0("********************* Filtered one. Only geneID ", rownames(queryData), " is in common."))
		return(res)
	}
	if (sd(queryData[,2]) == 0) {
		res <- getDefaultErrorDataframe("Uniform expression for common genes in signature. No correlation!") # no overlapping genes between query and reference
		print("********************* Uniform expression in signature for common genes. No correlation!")
		return(res)
	}
# 	if (dim(queryData)[1] == 1 | dim(queryData)[1] == 2) { ## Mehdi: it was only == 0 before!
# 		res$Remark <- paste0("Filtered few; not good for concordance. Only geneID(s) ", paste(rownames(queryData),collapse=","), " in common.") # no overlapping genes between query and reference
# 		print("********************* is it coming here inside?")
# 		return(res)
# 	}

	print("********************* is it coming here?")
#	refData[which(refData =="NaN")]<-NA
	# For Encode data eun GRS instead of CC 
# 	if(libName=="LIB_3" & dim(queryProfile)[2] == 3){
	if(libName=="LIB_3"){
	    if(dim(queryProfile)[2] == 3) {
#		source("/opt/raid10/genomics/software/RPrograms/source/batchTregGRS.R")
		
		##@Mehdi There is no need to source this file! 
		##        because "scaleGRS" function is not called neither here nor in "batchTregGRS" function
		
#		source("/opt/raid10/genomics/software/RPrograms/source/scaleGRS.R")

		#compute Up and Down Pvalues
		queryDataDown <- ifelse(queryData[,2]<0,queryData[,3]/2, 1-queryData[,3]/2)
		queryDataDown[queryDataDown<1e-10] <- 1e-10
		queryDataDown <- data.frame(queryData[,1],queryDataDown)

		queryDataUp <- ifelse(queryData[,2]>0,queryData[,3]/2, 1-queryData[,3]/2)
		queryDataUp[queryDataUp<1e-10] <- 1e-10
		queryDataUp <- data.frame(queryData[,1],queryDataUp)
		
		# Run TregGRS
		 grsZAll <- parallel::mclapply(1:dim(refData)[2],function(j) {
				print(j)
 				combineData <- na.omit(data.frame(EntrezID=queryDataUp[,1],Score=refData[,j],Prob=refProbsData[,j],queryDataUp[,2],queryDataDown[,2],stringsAsFactors=F))
				grsZUp <- batchTregGRS(query.p=combineData[,c(1,4)],reference.p=combineData[,c(1,2,3)] , na.rm = TRUE, estimateNullDistr = F, nullDistrQuantiles = c(0.9, 0.95, 0.99),
                               		nullDistrN = 100, tolerateWarnings=TRUE, pAdjust.method.query=NULL, pAdjust.method.reference=NULL, lambda=0.005,plotRescaled=F)[1:2]
     				grsZDown <- batchTregGRS(query.p=combineData[,c(1,5)],reference.p=combineData[,c(1,2,3)] , na.rm = TRUE, estimateNullDistr = F, nullDistrQuantiles = c(0.9, 0.95, 0.99),
                               		nullDistrN = 100, tolerateWarnings=TRUE, pAdjust.method.query=NULL, pAdjust.method.reference=NULL, lambda=0.005,plotRescaled=F)[1:2]
  #                 		grsZ <- ifelse( abs(grsZUp)> abs(grsZDown),grsZUp, grsZDown*(-1))
                   		grsP <- ifelse( abs(grsZUp[[1]])> abs(grsZDown[[1]]), grsZUp[[1]], ifelse(grsZDown[[1]] < 0, grsZDown[[1]], grsZDown[[1]]*(-1)))
						grsZ <- ifelse( abs(grsZUp[[2]])> abs(grsZDown[[2]]), grsZUp[[2]], ifelse(grsZDown[[2]] < 0, grsZDown[[2]], grsZDown[[2]]*(-1)))
   				list(grsZ, grsP)
 			}, mc.cores=8)
		grsZAll <- as.data.frame(data.table::rbindlist(grsZAll))												   
		rownames(grsZAll) <- colnames(refData)
		ord <- order(abs(as.numeric(grsZAll[,1])),decreasing=TRUE)
		#ord <- order(unlist(grsZAll),decreasing=TRUE)
                grsZAll <- grsZAll[ord,]
		#return All significant or top 10
		sigComp <- which(abs(as.numeric(grsZAll[,1])) >= 6)
		if(length(sigComp) == 0) { sigComp <- c(1:min(10, nrow(grsZAll))) }
		nGenes <- t(!is.na(queryData[,1])) %*% (!is.na(refData))												  
		corTable <- data.frame(dataset=rownames(grsZAll)[sigComp],CC=grsZAll[sigComp,1], PV=grsZAll[sigComp,2], nGenes=nGenes[1,sigComp], stringsAsFactors=F)
	    } else{
		res <- getDefaultErrorDataframe("No P-values provided")
		return(res)
	    }
	} else {
		# Run pair-wise correlation against RData and query profile
#@!		AllCor <- cor(queryData[,2],refData,use="pairwise.complete.obs")
# 		source("/opt/raid10/genomics/mehdi/ilincs/gitHub/cor_with_pval.R", local=T)
#@!		options(scipen=-999, digits=3)
		if (dim(queryData)[1] == 2) {
		    tmpcor <- cor(queryData[,2], refData, use="pairwise.complete.obs")
		    AllCor <- data.frame(dataset=colnames(tmpcor), CC=tmpcor[1,], PV=NA, nGenes=2, stringsAsFactors=FALSE)
		} else {		
		    AllCor <- cor_with_pval(queryData[,2], refData, cortype="pearson")#use="pairwise.complete.obs")
		}
		colnames(AllCor) <- c("dataset", "CC", "PV", "nGenes")
#@!		ord <- order(abs(AllCor), decreasing=TRUE)
		ord <- order(abs(AllCor$CC), decreasing=TRUE)
		#ord <- order(AllCor,decreasing=TRUE)
#@!		AllCor <- AllCor[ord]
		AllCor <- AllCor[ord,]
		AllCorNames <- colnames(refData)[ord]

		# Return the CC  results back to java to display the results
#@!		sigComp <- which(abs(unlist(AllCor)) > 0.1)
		sigComp <- which(abs(unlist(AllCor$CC)) > 0.1)
		if(length(sigComp) == 0){ # no significant comparisons then return top 10
			sigComp <- c(1:10)
		}
#@!		corTable <- data.frame(dataset=AllCorNames[sigComp],CC=AllCor[sigComp], stringsAsFactors=F)
		corTable <- AllCor[sigComp,]
	}
	# treeview of significant profiles
#	dataClust <- refData[,sigComp]
#	dataClust[which(dataClust =="NaN")]<-NA
#	ro <- apply(dataClust,2,function(x) sum(is.na(x)))
#	co <- apply(dataClust,1,function(x) sum(is.na(x)))
#	missingRows <- which(ro > (dim(dataClust)[1])/2)
#	missingCols <- which(co > (dim(dataClust)[2])/2) 
#  	dataClust <- dataClust[-missingRows,-missingCols] 
#	if((dim(dataClust)[1]) > 1){	
#		colClust<- hclust(dist(data.matrix(t(dataClust))),method="complete")
#		rowClust<- hclust(dist(data.matrix((dataClust))),method="complete")
#		colDendro <- as.dendrogram(colClust)
#		rowDendro <- as.dendrogram(rowClust)
#		create_treeview_files(res[[1]],dataClust,ifCluster="both",isCentered="FALSE",colClust,rowClust,path_to_write)
#	}
	
	res[[(length(res)+1)]] <- paste(corTable$dataset, collapse="DELIM")
        colnames(res)[dim(res)[2]] <- "tableNames"
        res[[(length(res)+1)]] <- paste(corTable$CC, collapse="DELIM") 
        colnames(res)[dim(res)[2]] <- "cc"
        res[[(length(res)+1)]] <- paste(corTable$PV, collapse="DELIM") 
        colnames(res)[dim(res)[2]] <- "pv"
        res[[(length(res)+1)]] <- paste(corTable$nGenes, collapse="DELIM") 
        colnames(res)[dim(res)[2]] <- "nGenes"
	# save it as excel file and send the path
	write.table(corTable,file=paste(path_to_write,"CorTable_",as.character(res[[1]]),".xls",sep=""), quote=TRUE, sep="\t", row.names=F, col.names=T, append=T)
        res[[(length(res)+1)]] <- paste(path_to_write,"CorTable_",as.character(res[[1]]),".xls",sep="") 
        colnames(res)[dim(res)[2]] <- "corTablePath"
        res$NoOfGenes <- length(genes)
#         if (is.na(res$NoOfProbes)) res$NoOfProbes <- res$NoOfGenes
        ## Mehdi: adding signatures meta data
        metacols <- switch(libName, 
			    LIB_1="datasetID,factor,Level1,Level2",
			    LIB_2="compound,concentration,cellLine",
			    LIB_3="treatment,antibodyTarget,cellLine",
			    LIB_5="compound,concentration,cellLine,time",
			    LIB_6="treatment,perturbagenID,cellLine,time",
			    LIB_7="treatment,perturbagenID,cellLine,time",
			    LIB_8="datasetID,factor,Level1,Level2",
			    LIB_10="perturbagenID,compound,cellLine")
#         sigmeta <- getSignatureMeta(prop=metacols, signatures=paste(as.vector(corTable[,1]),collapse=","), debuging=debuging, sigDB=sigDB)
	sigmeta <- getSignatureMeta(prop=metacols, signatures=paste(corTable$dataset, collapse=","), debuging=debuging, sigDB=sigDB)
	sigmatch <- match(corTable$dataset, sigmeta$signatureID)
	sigmeta <- sigmeta[sigmatch,]
	sigmeta <- sigmeta[,-1]
	if (libName=="LIB_1" | libName=="LIB_8") {
# 	    sigScores <- cbind(sigScores, sigmeta[,1:4])
# 	    colnames(sigScores) <- c("tableNames", "zScores", "datasetID","factor","Level1","Level2")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
 	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
 	colnames(res)[(dim(res)[2]-3):dim(res)[2]] <- c("datasetID","factor","Level1","Level2")
 	} else if  (libName=="LIB_2") {
# 	    sigScores <- cbind(sigScores, sigmeta[,1:3])
# 	    colnames(sigScores) <- c("tableNames", "zScores", "compound","concentration","cellLine")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
 	colnames(res)[(dim(res)[2]-2):dim(res)[2]] <- c("compound","concentration","cellLine")
 	} else if (libName=="LIB_3") {
# 	    sigScores <- cbind(sigScores, sigmeta[,1:3])
# 	    colnames(sigScores) <- c("tableNames", "zScores", "treatment","antibodyTarget","cellLine")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
 	colnames(res)[(dim(res)[2]-2):dim(res)[2]] <- c("treatment","antibodyTarget","cellLine")
 	} else if (libName=="LIB_5") {
# 	    sigScores <- cbind(sigScores, sigmeta[,1:4])
# 	    colnames(sigScores) <- c("tableNames", "zScores", "compound","concentration","cellLine","time")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
 	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
 	colnames(res)[(dim(res)[2]-3):dim(res)[2]] <- c("compound","concentration","cellLine","time")
	} else if (libName=="LIB_10") {
# 	    sigScores <- cbind(sigScores, sigmeta[,1:4])
# 	    colnames(sigScores) <- c("tableNames", "zScores", "compound","concentration","cellLine","time")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
 	colnames(res)[(dim(res)[2]-3):dim(res)[2]] <- c("perturbagenID","compound","cellLine")
 	} else {#	if (libName="LIB_6" | libName="LIB_7") {
# 	    sigScores <- cbind(sigScores, sigmeta[,1:4])
# 	    colnames(sigScores) <- c("tableNames", "zScores", "treatment","perturbagenID","cellLine","time")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,1]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,2]),collapse="DELIM")
	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,3]),collapse="DELIM")
 	res[[(length(res)+1)]] <- paste(as.vector(sigmeta[,4]),collapse="DELIM")
 	colnames(res)[(dim(res)[2]-3):dim(res)[2]] <- c("treatment","perturbagenID","cellLine","time")
 	}

# 	res <- list(usedGenes=usedGenes , sessionID=sessionID, Remark=remark, enrichTablePath=enrichTablePath, sigScores=sigScores)
 	#### Parsing 'res' for Volcano plot #####

# # # # #     if (reAnalyse) {
# # # # #       file.copy(paste(path_to_write,"temp_volcano_", previousID, ".xls", sep = ""), 
# # # # # 		paste(path_to_write,"temp_volcano_", as.character(res[[1]]), ".xls", sep = ""), overwrite = TRUE)
# # # # #       file.copy(paste(path_to_write,"temp_volcano_", previousID, ".RData", sep = ""), 
# # # # # 		paste(path_to_write,"temp_volcano_", as.character(res[[1]]), ".RData", sep = ""), overwrite = TRUE)
# # # # #       
# # # # #       } else parse_res(res, path_to_write, write=write)
# # # # #       parse_res(res, path_to_write, write=write)
      
        #### End parsing 'res' for Volcano plot #####
        
     	return(res)

     	} else { res <- createSignatureByDiffGEA(exp=exp, prop=prop,	filterchk=filterchk,	basement=basement,
					 includeORexclude=includeORexclude, ifCluster=ifCluster, up=up,	down=down, window.size=window.size, 
					 path_to_write=path_to_write,	pORqVal=pORqVal, userUploadedProfilePath=userUploadedProfilePath, previousID= previousID, reAnalyse=reAnalyse, 
					 glist=queryGenelist, pValue=pValue, foldchange=foldchange, MLinfo=MLinfo, write=write, homoloGenes=homoloGenes, multigroup=multigroup, 
					 sigDB=sigDB, chost=chost, debuging=debuging)

     	return(res)
     	}
}
