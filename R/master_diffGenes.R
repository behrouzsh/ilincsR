
#' A function for finding differentially expressed genes in a next-gene experiment.
#'
#' This function allows you to run differntial expression analysis in GenomicsPortals.
#' @param exp This argument is the experiment name, defined as a string corresponding to GenomicsPortals experiment name. The default is set to NULL.
#' @param prop prop defines which property of the experiment is going to be used for statistical analysis to specify the groups of samples.
#' @param filterchk filterchk is in the form of a string. It can have multiple property of the pData to filter
#'	based on. Each pair of property:value should be separated by a comma ",". Pairs are saparated
#'	by colon ":". An example ot filterchk is :
#'	filterchk="property1:value1,property1:value2,property2:value1"
#' @param basement A two value vector containing "baseline" and "treatment".
#' @param includeORexclude This argument basically is designed to filter the pData based on what 
#'	filterchk is or selected property of eset when is set to "1", or the reverse selection of filterchk 
#'	when set to "2". 
#'	It also can have NULL or "n" values when there is no property selected, in this case 
#'	the function will return the original ExpressionSet.
#' @param ifCluster This parameter defines in which way the data should be clustered. The default value in this work flow is set to "both" which 
#'	clusters both genes (rows) and samples (columns) in the data matrix.
#' @param pValue pValue is a numeric value between 0 and 1 and determines the size of the tests (ttest/ftest or the alpha value) 
#'	in selecting differentially expressed genes.
#' @param foldchange This is also a numeric value for selecting differentially expressed genes based on their log2 fold changes between groups.
#' 	If set to NULL along with pValue of NULL the function will use top 100 genes based on p-values. 
#' @param up,down Specifies the maximum number of up and down regulated genes.
#' @param window.size ...
#' @param path_to_write This necessary parameter specifies the path which user wants to save the 
#'	output files (heatmaps, TreeView files and some other output files) to be saved in.
#' @param pORqVal is a string argument and defines if "p" values are used in the analysis or "q" values.
#' @param previousID This will be used if user wants to reanalyse a new dataset or experiment using the gene list derived from first round of analysis.
#' @param reAnalyse Should be set to TRUE only if same/another experiment or dataset is needed to be analysed based on the results of the first round of analysis.
#' @param MLinfo This argument is a list object containing previous analysis arguments for the function and should be used only if 
#'	reAnalyse is set to TRUE. After each round of running "master_LincsDataAnalysis" function, all the function arguments are saved
#'	just in case user wants to proceed any reAnalysing with a new dataset.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param multigroup A critical parameter for specifying which work flow user wants to choose. The default is set to FALSE meaning the function
#'	will run in the two group mode and will calculates the signatures in specified groups (only two groups are allowed). 
#'	In the other hand, if multigroup is set to TRUE, it will run an F test between all selected groups (selected previously in filterchk).
#' @param authors M. Fazel-Najafabadi
#' @keywords exp, prop, reAnalyse.
#' @export 
#' @examples
#' ## not run
#' res_multigroup <- master_diffGenes(exp="EDS-1014", prop="subtype", filterchk="subtype:Basal,subtype:Luminal,subtype:Claudin-low", 
#'	basement=NULL, includeORexclude="1", ifCluster="both", pORqVal="p", pValue=0.000001, foldchange=4, 
#'	up=4000, down=1000, window.size=50, path_to_write="/patho/to/your/output/directory/", multigroup=TRUE)
#'
#' # For reanalysing:
#' res2_multigroup <- master_diffGenes(exp="EDS-1014", prop="subtype", filterchk="subtype:Basal,subtype:Luminal,subtype:Claudin-low", 
#' basement=NULL, includeORexclude="1", ifCluster="both", pORqVal="p", pValue=0.000001, foldchange=4, 
#' up=4000, down=1000, window.size=50, path_to_write="/patho/to/your/output/directory/", multigroup=TRUE, 
#' previousID="find previousID from your saved results in the first step", reAnalyse=TRUE, MLinfo=MLinfo)
#'
#'	# previousID may look something like this: "Mon_Oct_26_23_49_56_2015_1084923"
#'	# "Date + a randomely generated number separated by under scores"
#' ## not run
#' @export 

master_diffGenes <- function(exp, prop, filterchk="n", basement=NULL, includeORexclude, ifCluster="both", up=4000, down=1000, window.size=50, 
			    path_to_write, pORqVal, previousID=NULL, reAnalyse=FALSE, glist=NULL, pValue=NULL, foldchange=NULL, 
			    MLinfo=NULL, homoloGenes=TRUE, multigroup=FALSE, sigDB="ilincs_sigs", chost="ilincs.org", debuging=FALSE) 
{

if(is.null(filterchk)) filterchk <- "n"
if(is.null(prop)) prop <- "n"
filterchk <- gsub("<NA>", "NA", filterchk)
# # # # # if (!multigroup) {
	## log
	a <- Sys.time()
	print(paste("master_diffGenes called with following parameters at",a))
        print(paste("previousID",previousID))
        print(paste("exp",exp))
        print(paste("prop",prop))
        print(paste("filterchk",filterchk))
        print(paste("basement",basement))
        print(paste("includeORexclude",includeORexclude))
       # print(paste("pValueCutoff",pValue))
       # print(paste("foldchange",foldchange))
        print(paste("pORqVal",pORqVal))
        print(paste("path_to_write",path_to_write))
        print(paste("reAnalyse",reAnalyse))
#	print(paste("before get pdata"))
	print(paste("before get pdata changed"))

	startQuery <- Sys.time()
##	library(RMySQL)
##	source("http://eh3.uc.edu/r/GRS_prod.R")
##	source("http://eh3.uc.edu/r/mod_Rserve_prod.R")
#	source("/opt/raid10/genomics/mehdi/ilincs/code/R/iLincs_R_Scripts/GRS_prod.R")
#	source("/opt/raid10/genomics/mehdi/ilincs/code/R/iLincs_R_Scripts/mod_Rserve_prod.R")

#	source("/opt/raid10/genomics/software/RPrograms/source/ilincs/dev/GRS_prod.R") # new locations
#	source("/opt/raid10/genomics/software/RPrograms/source/ilincs/dev/mod_Rserve_prod.R")	

	## Get dataType of the experiment
# # # # # 	servSet <- getServerSettings(test)
# # # # # 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#@!	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
        
#        if (!test) {
#        mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="10.165.4.231", port = 4040, password="public")
#        } else {
#        mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="gimm2.ketl.uc.edu", password="public")}
	print("db changed")
        #mycon <- dbConnect(MySQL(), user="public", dbname="QueryTables", host="db2.ketl.uc.edu", port = 3306, password="public")
# # # # #         sql <- paste("select Platform,DataType,DataFormat,eset_path from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#@!        sql <- paste("select platformDb,DataType,DataFormat,eset_path from ExperimentMetadata where experimentName = '",exp,"'",sep="")
# # # # #         rs <- DBI::dbSendQuery(mycon, sql)
# # # # #         dt <- DBI::fetch(rs, n = -1)
# # # # #         dataType <- dt[1,"DataType"]
# # # # #         dataFormat <- dt[1,"DataFormat"]
# # # # #         platform <- dt[1,"Platform"]
# # # # #         esetpath <- dt[1,"eset_path"]
# # # # #         esetpath <- try(sub("mnt", "opt", esetpath))
# # # # #         DBI::dbDisconnect(mycon)
        info <- expInfo(exp, debuging)
        dataType <- info["DataType"]
        dataFormat <- info["DataFormat"]
        platform <- info["Platform"]
        esetpath <- info["eset_path"]
        esetpath <- try(sub("mnt", "opt", esetpath))
	db <- platform
	## Get pData using getDataForGRS function in GRS.R
# # # # # 	pData <- getpData(exp, platform, test)
# # # # # 	rownames(pData) <- gsub(paste0(exp, "_"), "", pData$MeasurementName) ## Mehdi: Added 2018 for matching samples
# # # # # 	pData[is.na(pData)] <- "NA"
	
	print("after get pdata")

	## Get expression data using getDataForGRS function in GRS.R
       	esetName <- paste(exp,"_",platform,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
       	## Mehdi: added the print
       	print(paste(esetName))

	## Mehdi: New option for gene pipeline
	if(!is.null(pValue)) pValue <- as.numeric(pValue)
	if(!is.null(glist) & is.null(pValue)) {
		    eset <- try(getData(exp=exp, glist=glist, filterchk=filterchk, includeORexclude=includeORexclude, homoloGenes=homoloGenes, debuging=debuging))
# 		    if(length(grep("^Error.*",eset))){  # if eset not found at the specified location
		    if(class(eset) != "ExpressionSet"){  # if eset not found at the specified location
			    esetpath=""  # query from database
		    } else {
# 			    eset <- as(get(eset),"ExpressionSet")
# 			    if(sum(is.na(match(rownames(pData),colnames(Biobase::exprs(eset))))) > 0)  # not matching
# 			    {
# 				    colnames(Biobase::exprs(eset)) <- rownames(pData)
# 			    }
# 			    Biobase::pData(eset) <- pData[colnames(Biobase::exprs(eset)), , drop=FALSE]
# 			    fdat <- annotateFeatures(exp, debuging=debuging)
# 			    fdat <- fdat[match(rownames(Biobase::exprs(eset)), fdat$PROBE), ]
# 			    rownames(fdat) <- fdat$PROBE
# 			    Biobase::fData(eset) <- fdat
			    assign(esetName,eset)
			    probes <- rownames(Biobase::exprs(eset))
			    geneprobes <- get_gene_probe_eset(db,unique(probes), debuging)
		    }
	} else {
	    if(esetpath != "")  # esetpath present. directly load the eset from the directory
	    {
		    eset <- try(load(esetpath))
		    if(length(grep("^Error.*",eset))){  # if eset not found at the specified location
			    esetpath=""  # query from database
		    } else {
			    eset <- as(get(eset),"ExpressionSet")
# # # # # 			    if(sum(is.na(match(rownames(pData),colnames(Biobase::exprs(eset))))) > 0)  # not matching
# # # # # 			    {
# # # # # # 				    colnames(Biobase::exprs(eset)) <- rownames(pData)
# # # # # 				Biobase::exprs(eset) <- Biobase::exprs(eset)[ , rownames(pData), drop=FALSE]
# # # # # 			    }
# # # # # 			    Biobase::pData(eset) <- pData[colnames(Biobase::exprs(eset)), , drop=FALSE]
			    fdat <- annotateFeatures(exp, debuging=debuging)
# 			    fdat <- fdat[match(rownames(Biobase::exprs(eset)), fdat$PROBE), ]
			    fdat <- fdat[rownames(Biobase::exprs(eset)), ]
# 			    rownames(fdat) <- fdat$PROBE
			    Biobase::fData(eset) <- fdat
			    assign(esetName,eset)
			    probes <- rownames(Biobase::exprs(eset))
			    geneprobes <- get_gene_probe_eset(db,unique(probes), debuging)
		    }
#       	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
	    
	    }
	    if(esetpath == "")  # eset path info not provided then obtain the data by querying the database
	    {
		    eData<-getDataForGRS(exp,filterchk,basement,includeORexclude,samples=rownames(pData),platform)
		    eData<-eData[which(eData[,"Gene"] != "" ),]

		    ## Create eset
#			ord <- order(pData[,prop])
#       		pData <- pData[ord,]
#        		NoOfProbes <- dim(eData)[1]/dim(pData)[1]
#       		geneNames <- eData[1:NoOfProbes,"Gene"]

		    ExprData <- NULL
		    ExprData<-matrix(eData[,"Value"],ncol=dim(pData)[1])
		    # unqVals <- as.character(unique(eData[,1]))
		    # unqIndex <- match(unqVals,rownames(pData))

		    rownames(ExprData) <-eData[1:length(unique(eData[,"Probe"])),"Probe"] 
		    colnames(ExprData)<-rownames(pData)
		    ExprData<-as.data.frame(ExprData)

		    ## Call functions from plot_with_properties_eset wrapper from mod_Rserve.R
		    probes<-rownames(ExprData)
		    geneprobes<-get_gene_probe_eset(db,unique(probes), debuging)
		    ExprData<-ExprData[geneprobes[,"probe_name"],]
		    
		    ## need to add fData later if needed for this pipeline!

##		    library(Biobase)
		    assign(esetName,new("ExpressionSet", exprs=as.matrix(ExprData),
					phenoData=new("AnnotatedDataFrame",data=pData),
					featureData=new("AnnotatedDataFrame",data=as.data.frame(rownames(ExprData), row.names=rownames(ExprData), stringsAsFactors=FALSE))))

		    endQuery <- Sys.time()
		    print(paste("start time for query and eset",startQuery))
		    print(paste("end time for query and eset",endQuery))
		    print(paste("total time for query and eset",endQuery - startQuery))
#       	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
	    } # if no esetpath found
	}

	## filter eset according to criterion and then run anova
	    f <- filter_eset(get(esetName),filterchk,basement,includeORexclude)
#	    save(f, file=paste(path_to_write, "/myfiltered_eset.RData", sep="")) # Mehdi: for NA testing
	    
#  	remove NA  or missing data
	    if(!is.null(glist)) geneprobes <- geneprobes[geneprobes$gene %in% unlist(strsplit(glist, ",")),]
	    ExprData <- Biobase::exprs(f)[rownames(geneprobes), , drop=FALSE]
	    countingMissingData <- apply(ExprData,1,function(x) sum(is.na(x)))
	    nonMissingData <- (countingMissingData != length(colnames(ExprData)))
	    ExprData <- ExprData[nonMissingData, ,drop=FALSE]
	    ExprData <- as.matrix(ExprData)
# 	    rownames(ExprData) <- as.vector(rownames(ExprData)[names(nonMissingData)])

	    pData <- Biobase::pData(f)
	    pData[is.na(pData)] <- "NA"
	    if(!is.null(basement)) {
		baseline <- pData[,prop] == basement[1]
		treatment <- !baseline
		pData <- rbind.data.frame(pData[baseline,,drop=F], pData[treatment,,drop=F], stringsAsFactors = FALSE)
		ExprData <- ExprData[,rownames(pData), drop=F]
	    }

	    fData <- Biobase::fData(f)[rownames(ExprData), , drop=FALSE]
	    fData <- as.data.frame(apply(fData, 2, as.character), stringsAsFactors=FALSE, row.names=rownames(fData))
	    ## Set prop for "n"
	    if(prop == "n") {
		    pr <- apply(pData, 2, function(j) length(unique(j)))
		    prp <- colnames(pData)
		    prp <- prp[!pr==1]
		    pr <- pr[!pr==1]
		    prop <- prp[pr==min(pr)][1]
	    }
	    ## Mehdi for reAnalysing old method
# 	    if (reAnalyse) {
# 		print(paste("...reAnalysing"))
# # 		fromVolcano <- read.table(paste(path_to_write, "genelist_temp_volcano_", previousID, ".xls", sep=""), sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
# 		fromVolcano <- get(load(paste0(path_to_write, "temp_volcano_", previousID, ".RData")))
# 		colnames(fromVolcano) <- c("geneID","geneName","coefficients","Pvals") # c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
# 		if(!is.null(pValue)) {
# 		    fromVolcano <- fromVolcano[(!is.na(fromVolcano$Pvals) & fromVolcano$Pvals <= as.numeric(pValue)), ]
# 		}
# 		if(!is.null(foldchange)) {
# 		    foldchange <- as.numeric(unlist(strsplit(as.character(foldchange), "DELIM")))
# 		    if(length(foldchange)==1) fromVolcano <- fromVolcano[(fromVolcano$coefficients <= -abs(foldchange) | fromVolcano$coefficients >= abs(foldchange)), ]
# 		    if(length(foldchange)==2) fromVolcano <- fromVolcano[(fromVolcano$coefficients <= -abs(foldchange[1]) | fromVolcano$coefficients >= abs(foldchange[2])), ]
# 		}
# 		geneprobes <- geneprobes[match(fromVolcano$geneID, geneprobes$gene),]
# 		geneprobes <- unique(geneprobes)
# # 		geneprobes <- geneprobes[geneprobes$gene %in% fromVolcano$GeneID ,]
# 		ExprData <- Biobase::exprs(f)[geneprobes[!is.na(geneprobes$probe_name),"probe_name"], ] #Biobase::exprs(f)[geneprobes[,"probe_name"],]
# # 		ExprData <- ExprData[unique(rownames(ExprData)), ]
# 	    }
       	## Mehdi: for reAnalysing
#       	} else {
#	    load(paste(path_to_write, "filteredesetinfo_", previousID, ".RData", sep=""))
#	    f <- filter_eset(saved.eset, filterchk,basement,includeORexclude)
#	    ExprData <- Biobase::exprs(f)
#	    pData <- Biobase::pData(f)
#	}
	    
	# check for 2 levels for comparison
	
	if(!is.null(basement)) {
	    status <- factor(as.character(pData[,prop]), levels=basement)
	} else {
	    status <- as.factor(as.character(pData[,prop]))
	}
# # # # # 	if(length(levels(status)) != 2){
# # # # #                 remark <- "No two levels"
# # # # #                 gpgene <- "NA"
# # # # #                 gpprobe <- "NA"
# # # # #                 dc2 <- "NA"
# # # # #                 pkeggvalue <- "NA"
# # # # #                 sample_after_filtering <- "NA"
# # # # #                 filteredSamples <- "NA"
# # # # #                 sessionID <- "NA"
# # # # #                 pvalueLR <- "NA"
# # # # #                 noOfGenes <- "NA"
# # # # #                 noOfProbes <- "NA"
# # # # #                 coeffs <- "NA"
# # # # #                 coeffsGeneID <- "NA"
# # # # # 		coeffsGeneName <- "NA"
# # # # #                 sample_after_filtering <- 0
# # # # #                 filteredSamples <- "NA"
# # # # #                 level1 <- paste(as.vector(levels(status)),collapse="DELIM")
# # # # #                 level2 <- "NA"
# # # # # 		Pvals <- "NA"
# # # # #                 x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR,noOfGenes,noOfProbes,coeffs,coeffsGeneID,coeffsGeneName,level1,level2,Pvals)
# # # # #                 colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR","NoOfGenes","NoOfProbes","coefficients","coefficientsGeneID","coefficientsGeneName","level1","level2","Pvals")
# # # # #                 return(x)
# # # # # 	}

	print("before Anova")
	## Step 1 compute anova

	startanova <- Sys.time()
##	library(limma)
##	library(genefilter)
        #design <- model.matrix(~ -1+status)
        design <- model.matrix(~status)
        fit <- limma::lmFit(ExprData, design = design)
 #      junk<-summary(lm(ExprData[1,]~status))
 #	junk<-summary(lm(ExprData[1,]~-1+status))
        fit2 <- NULL
	try(fit2 <- genefilter::rowFtests(ExprData,status))
#       try(fit2 <- eBayes(fit),TRUE)

	if(is.null(fit2))
        {
                #return(paste("Error in anova",geterrmessage()))
		print("Error in anova")
		remark <- "Error in anova"
		gpgene <- "NA"
		gpprobe <- "NA"
		dc2 <- "NA"
		pkeggvalue <- "NA"
		sample_after_filtering <- "NA"
		filteredSamples <- "NA"
		sessionID <- "NA"
		pvalueLR <- "NA"
		noOfGenes <- "NA"
                noOfProbes <- "NA"
                coeffs <- "NA"
                coeffsGeneID <- "NA"
		coeffsGeneName <- "NA"
                sample_after_filtering <- 0
                filteredSamples<-"NA"
		level1 <- "NA"
		level2 <- "NA"
		Pvals <- "NA"
                x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR,noOfGenes,noOfProbes,coeffs,coeffsGeneID,coeffsGeneName,level1,level2,Pvals)
                colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR","NoOfGenes","NoOfProbes","coefficients","coefficientsGeneID","coefficientsGeneName","level1","level2","Pvals")
		print(paste("Current sessionID is: ", x$sessionID))

                return(x)
        }
        
	fit2 <- fit2[complete.cases(fit2), ,drop=FALSE]
	fit <- fit$coefficients[rownames(fit2), ,drop=FALSE]
        
        
	print("after Anova")
        ## if exactly 2 levels use p.values else use F.p.value
# # # # #         if(length(unique(pData[,prop])) == 2) {
# # # # #                 #Pvalues <- fit2$p.value[,2]
# # # # #                 Pvalues <- fit2$p.value
# # # # #                 names(Pvalues)<-rownames(fit2)
# # # # #         } else {
                #Pvalues <- fit2$F.p.value
                #names(Pvalues)<-names(fit2$p.value[,2])
                Pvalues <- fit2$p.value
                names(Pvalues)<-rownames(fit2)
# # # # #         }

	if(pORqVal=="q"){
##		library(qvalue)
		Pvalues <- qvalue::qvalue(p=Pvalues,lambda=0)$qvalues
	}

        if(length(unique(pData[,prop])) == 2) {
# # # # # 		coeffs <- fit$coefficients[,2]
		coeffs <- fit[,2]
		multigroup <- FALSE
	} else {  # if more than 2 levels take the mx fold change
# # # # # 		lt <- fit$coefficients
#		lt1 <- lt[match(rownames(ExprData),rownames(lt)),]
# # # # #         	coeffs <- apply(lt[,-1],1,function(x) max(marray::rm.na(x)))
        	coeffs <- apply(fit[,-1],1,function(x) max(marray::rm.na(x)))
        	multigroup <- TRUE
	}
	
	## order  data according to Pvalues
# # # # # 	ord <- order(Pvalues)
# # # # # 	Pvalues <- Pvalues[ord]
# # # # # 	ExprData <- ExprData[ord,]
# # # # # 	coeffs <- coeffs[ord]
# # # # # 	fit <- fit[ord, ]
	Pvalues <- Pvalues[order(Pvalues)]
	Pvalues <- Pvalues[!is.na(Pvalues)]
	ExprData <- ExprData[match(names(Pvalues),rownames(ExprData)),]
	coeffs <- coeffs[match(names(Pvalues),names(coeffs))]
# 	fit <- fit[names(Pvalues), ]
	geneprobes <- geneprobes[match(names(Pvalues),rownames(geneprobes)), ]
	completeSig <- geneprobes[, c(2,1,3,4)]
	colnames(completeSig) <- c("PROBE", "ID_geneid", "Name_GeneSymbol", "Description")
	rownames(completeSig) <- completeSig$PROBE

	perPropertyAverages <- NULL
	factorLevels <- unique(pData[,prop])
# # # # # 	if(filterchk != "n") {
# # # # # 	    flch <- unlist(strsplit(filterchk, ",,,"))
# # # # # 	    flch <- flch[grep(prop, flch)]
# # # # # 	    flch <- gsub(prop, "", flch)
# # # # # 	    flch <- gsub(":", "", flch)
# # # # # 	    factorLevels <- flch
# # # # # 	}
	for(i in 1:length(factorLevels)) {
		temp <- which(as.vector(colnames(ExprData)) %in% rownames(pData)[pData[,prop]==factorLevels[i]])
		perPropertyAverages <- cbind(perPropertyAverages, apply(as.matrix(ExprData[,temp]),1,function(x) mean(x,na.rm=TRUE)))
	}
	colnames(perPropertyAverages) <- factorLevels
	completeSig <- cbind(completeSig, perPropertyAverages)

	completeSig$Value_LogDiffExp <- coeffs
	completeSig$Significance_pvalue <- Pvalues

# # # # # 	completeSig <- completeSig
# # # # # 	ExprData1 <- ExprData
# # # # # 	Pvalues1 <- Pvalues
# # # # # 	coeffs1 <- coeffs

# # # # # 	    if (reAnalyse) {
# # # # # 		print(paste("...reAnalysing"))
		if(!is.null(pValue)) {
		    Pvalues <- Pvalues[!is.na(Pvalues)]
		    Pvalues <- Pvalues[Pvalues <= pValue]
##		    Pvalues <- Pvalues[order(Pvalues)]  ## extra ordering
		    coeffs <- coeffs[match(names(Pvalues),names(coeffs))]
		    ExprData <- ExprData[match(names(Pvalues),rownames(ExprData)), , drop=FALSE]
		}
		if(!is.null(foldchange)) {
		    foldchange <- as.numeric(unlist(strsplit(as.character(foldchange), "DELIM")))
		    if(length(foldchange)==1) coeffs <- coeffs[(coeffs <= -abs(foldchange) | coeffs >= abs(foldchange))]
		    if(length(foldchange)==2) coeffs <- coeffs[(coeffs <= -abs(foldchange[1]) | coeffs >= abs(foldchange[2]))]
		    Pvalues <- Pvalues[match(names(coeffs),names(Pvalues))]
		    Pvalues <- Pvalues[order(Pvalues)]  ## extra ordering
		    coeffs <- coeffs[match(names(Pvalues),names(coeffs))]	## extra ordering
		    ExprData <- ExprData[match(names(coeffs),rownames(ExprData)), , drop=FALSE]
		}
# 		geneprobes <- geneprobes[match(fromVolcano$geneID, geneprobes$gene),]
# 		geneprobes <- unique(geneprobes)
# # 		geneprobes <- geneprobes[geneprobes$gene %in% fromVolcano$GeneID ,]
# 		ExprData <- Biobase::exprs(f)[geneprobes[!is.na(geneprobes$probe_name),"probe_name"], ] #Biobase::exprs(f)[geneprobes[,"probe_name"],]
# # 		ExprData <- ExprData[unique(rownames(ExprData)), ]
# # # # # 	    }
	
	endanova <- Sys.time()
	print(paste("start time for anova",startanova))
	print(paste("end time for anova",endanova))
	print(paste("total time for anova",endanova - startanova))
	geneprobes <- geneprobes[complete.cases(geneprobes[names(Pvalues),]),]
# # # # # 	ExprData2 <- ExprData
# # # # # 	Pvalues2 <- Pvalues
# # # # # 	coeffs2 <- coeffs

	if(dim(ExprData)[1] > 0)
	{
	# Get genes
# # # # # 	x <- match(rownames(ExprData),geneprobes[,"probe_name"])
# # # # # 	ExprData <- ExprData[which(is.na(x) == FALSE),]
# # # # # 	Pvalues <- Pvalues[which(is.na(x) == FALSE)]
# # # # # 	coeffs <- coeffs[which(is.na(x) == FALSE)]
# # # # #         geneprobes <- geneprobes[x[which(!is.na(x))],]	

#	# Get genes
#	geneprobes<-geneprobes[match(names(Pvalues),geneprobes[,"probe_name"]),]
#	x<-match(geneprobes[,"probe_name"],rownames(ExprData))
#	if(sum(is.na(x)) > 0)
#	{
#		x <- x[which(is.na(x) == FALSE)]
#	#	colnames(ExprData) <- paste(exp,"_",colnames(ExprData),sep="" )
#
#	}
#	ExprData<-ExprData[x,]

	# Create eset with only first 100 genes
        esetName <- paste(exp,"_",platform,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
        esetName <- gsub("-", "_", esetName)
##        assign(esetName, new("ExpressionSet",exprs=as.matrix(ExprData[1:100,]),phenoData=new("AnnotatedDataFrame",data=pData)))
	## Mehdi : eset for save
	esetName4save <- paste(esetName, "4save", sep="")
##        assign(esetName4save, new("ExpressionSet",exprs=as.matrix(ExprData),phenoData=new("AnnotatedDataFrame",data=pData)))
	
#	# pass the computed Pvalues
#	Pvalues <- Pvalues[match(rownames(ExprData),names(Pvalues))]
#	coeffs <- coeffs[match(rownames(ExprData),names(coeffs))]

	# Do clustering of only first 100 genes 
	ifCluster <- "both"
	# LR not applicable
	ifLR <- 0
	startA <- Sys.time()
# # # # # 	fData <- fData[match(rownames(ExprData), rownames(fData)), ]
	fData <- fData[rownames(ExprData), ]
	## Mehdi added for reanalyse  as.data.frame(rownames(ExprData), row.names=rownames(ExprData), stringsAsFactors=FALSE)
# # # # # 	if (!reAnalyse) {
	    if(is.null(glist) & is.null(pValue)) {
	      norows <-  1:min(100, nrow(ExprData)) 
	    } else {
	      norows <-  1:nrow(ExprData)
	    }
# # # # # 	    if(is.null(glist) & is.null(pValue)) {
# # # # # 	      norows <- names(Pvalues)[1:min(100, nrow(ExprData))] # 1:min(100, nrow(ExprData)) 
# # # # # 	    } else {
# # # # # 	      norows <- names(Pvalues)[1:nrow(ExprData)] # 1:nrow(ExprData)
# # # # # 	    }
            
# 	subsetSig <- completeSig[na.omit(match(rownames(ExprData)[norows],rownames(completeSig))), ]
	subsetSig <- completeSig[rownames(ExprData)[norows], ]
# 	subsetSig <- subsetSig
	
            assign(esetName, new("ExpressionSet", exprs=as.matrix(ExprData[norows,]), 
				  phenoData=new("AnnotatedDataFrame",data=pData), 
				  featureData=new("AnnotatedDataFrame",
# 					data=as.data.frame(rownames(ExprData[1:min(100, nrow(ExprData)),]), row.names=rownames(ExprData[1:min(100, nrow(ExprData)),]), stringsAsFactors=FALSE))
					data=fData[norows, ])
				  ))
# # # # #             assign(esetName4save, new("ExpressionSet", exprs=as.matrix(ExprData), 
# # # # # 				  phenoData=new("AnnotatedDataFrame",data=pData),
# # # # # 				  featureData=new("AnnotatedDataFrame",
# # # # # # 					data=as.data.frame(rownames(ExprData[1:min(100, nrow(ExprData)),]), row.names=rownames(ExprData[1:min(100, nrow(ExprData)),]), stringsAsFactors=FALSE))
# # # # # 					data=fData)
# # # # # 				  ))

# # # # #             if(is.null(glist) & is.null(pValue)) {
		
#           resultdf <- try(plot_with_properties_eset(eset=get(esetName), geneprobes[norows,], prop, web=FALSE, dataType,
# 						filterchk, basement, includeORexclude, ifCluster, ifLR, exp, db=platform, path_to_write,
# 						dataFormat, up, down, window.size, diffGenesSelected=TRUE, Pvalues[norows],
# 						eset4save=get(esetName4save), allPvalues=Pvalues, allcoeffs=coeffs,
# 						allgeneprobes=geneprobes, MLinfo=MLinfo))
		resultdf <- try(plot_with_properties_eset(eset=get(esetName), geneprobes[rownames(ExprData)[norows],], prop, web=FALSE, dataType, # geneprobes[rownames(ExprData)[norows],]
				filterchk, basement, includeORexclude, ifCluster, ifLR, exp, db=platform, path_to_write,
				dataFormat, up, down, window.size, diffGenesSelected=TRUE, Pvalues[norows],
				eset4save=get(esetName), allPvalues=NULL, allcoeffs=NULL,
				allgeneprobes=NULL, MLinfo=NULL, sigDB=sigDB, chost=chost, debuging=debuging))
# # # # #             } else {
# # # # # 		resultdf <- try(plot_with_properties_eset(eset=get(esetName4save), geneprobes[rownames(ExprData),], prop, web=FALSE, dataType,
# # # # # 				filterchk, basement, includeORexclude, ifCluster, ifLR, exp, db=platform, path_to_write,
# # # # # 				dataFormat, up, down, window.size, diffGenesSelected=TRUE, Pvalues,
# # # # # 				eset4save=get(esetName4save), allPvalues=NULL, allcoeffs=NULL,
# # # # # 				allgeneprobes=NULL, MLinfo=NULL))
# # # # #             }
	    endA <- Sys.time()
	    print(paste("start time for plot_with_properties_eset",startA))
	    print(paste("stop time for plot_with_properties_eset",endA))
	    print(paste("total time for plot_with_properties_eset",endA-startA))
	    print(paste("done",resultdf$Remark[1],"at time",Sys.time()))

	    resultdf[[(length(resultdf)+1)]] <- paste(length(unique(geneprobes[norows,"gene"])))
	    colnames(resultdf)[dim(resultdf)[2]] <- "NoOfGenes"
	    
	    resultdf[[(length(resultdf)+1)]] <- paste(length(geneprobes[norows,"gene"])) 
	    colnames(resultdf)[dim(resultdf)[2]] <- "NoOfProbes"
	    
	    resultdf[[(length(resultdf)+1)]] <- paste(as.vector(coeffs),collapse="DELIM")
	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficients"

	    resultdf[[(length(resultdf)+1)]] <- paste(geneprobes[match(names(coeffs),geneprobes[,"probe_name"]),"probe_name"],collapse="DELIM")
	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficientsProbe"
	    
	    resultdf[[(length(resultdf)+1)]] <- paste(geneprobes[match(names(coeffs),geneprobes[,"probe_name"]),"gene"],collapse="DELIM")
	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficientsGeneID"
	    
	    resultdf[[(length(resultdf)+1)]] <- paste(geneprobes[match(names(coeffs),geneprobes[,"probe_name"]),"symbol"],collapse="DELIM")
	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficientsGeneName"
	    
	if(!multigroup) {
	    
	    resultdf[[(length(resultdf)+1)]] <- levels(status)[1] 
	    colnames(resultdf)[dim(resultdf)[2]] <- "level2"
	    
	    resultdf[[(length(resultdf)+1)]] <- levels(status)[2] 
	    colnames(resultdf)[dim(resultdf)[2]] <- "level1"
	    
	    resultdf[[(length(resultdf)+1)]] <- length(status[status==resultdf[1,"level2"]]) 
	    colnames(resultdf)[dim(resultdf)[2]] <- "level2no"
	    
	    resultdf[[(length(resultdf)+1)]] <- length(status[status==resultdf[1,"level1"]]) 
	    colnames(resultdf)[dim(resultdf)[2]] <- "level1no"
	}
	    
	    resultdf[[(length(resultdf)+1)]] <- paste(as.vector(Pvalues),collapse="DELIM")
	    colnames(resultdf)[dim(resultdf)[2]] <- "Pvals"
	    
	    ## Mehdi: I added for volcano creation
	    #SignatureTable$coefficients <- coeffs
	    #SignatureTable$Pvals <- Pvalues
	    #.frame("GeneID"=dataClust[,1], ""=dataClust[,2], "Pvals"=Pvalues, stringsAsFactors = F)
	    #volcano_name = paste(path_to_write,"temp_volcano2_", as.character(sessionID), ".RData", sep = "")
	    #volcano_name2 = paste(path_to_write,"temp_volcano2_", as.character(sessionID), ".txt", sep = "")
	    #write.table(SignatureTable, file = volcano_name2, row.names=FALSE, sep="\t")
	    #save(SignatureTable, file = volcano_name)
	
	    print(paste("Current sessionID is: ", resultdf$sessionID))
# # # # # 	    write.table(completeSig, file=paste(path_to_write, "fullSig_", resultdf$sessionID, ".xls", sep=""), sep="\t", quote=FALSE, row.names=F, col.names=T, append=T)
# # # # # 	    write.table(subsetSig, file=paste(path_to_write, "subSig_", resultdf$sessionID, ".xls", sep=""), sep="\t", quote=FALSE, row.names=F, col.names=T, append=T)
	    write.table(completeSig, file=paste(path_to_write, "completeSig_", resultdf$sessionID, ".xls", sep=""), sep="\t", quote=TRUE, row.names=F, col.names=T, append=T)
	    SignatureTable <- completeSig
	    save(SignatureTable, file=paste(path_to_write, "temp_volcano_", resultdf$sessionID, ".RData", sep=""))
	    write.table(subsetSig, file=paste(path_to_write, "subsetSig_", resultdf$sessionID, ".xls", sep=""), sep="\t", quote=TRUE, row.names=F, col.names=T, append=T)
	    
	saved.eset <- get(esetName)
	means <- apply(data.matrix(Biobase::exprs(saved.eset)),1,median,na.rm=T)
	Biobase::exprs(saved.eset) <- as.matrix(sweep(data.matrix(Biobase::exprs(saved.eset)),1,means,"-"))
	Biobase::fData(saved.eset) <- Biobase::fData(saved.eset)[rownames(Biobase::exprs(saved.eset)), ]
	pd <- Biobase::pData(saved.eset)
	Biobase::pData(saved.eset) <- pd[colnames(Biobase::exprs(saved.eset)), ,drop=FALSE]
	rownames(Biobase::pData(saved.eset)) <- paste(Biobase::pData(saved.eset)[ , prop], rownames(Biobase::pData(saved.eset)), sep=".")
	save(saved.eset, file=paste(path_to_write, "/filteredeset_", resultdf$sessionID,".RData", sep=""))
	esetToGct1.3(saved.eset, paste0("/filteredeset_", resultdf$sessionID), path_to_write)
	    
	    
	    return(resultdf)
# # # # # 	} else {
# # # # # 
# # # # # 	    print(paste(path_to_write, "filteredesetinfo_", previousID, ".RData", sep=""))
# # # # # 	    print(paste(path_to_write, "genelist_temp_volcano_", previousID, ".xls", sep=""))
# # # # # 	    print(paste("********************* reAnalysing! "))
# # # # # 	    
# # # # # #	    load(paste(path_to_write, "filteredesetinfo_", previousID, ".RData", sep=""))
# # # # # #	    fromVolcano <- read.csv(paste(path_to_write, "genelist_temp_volcano_", previousID, ".csv", sep=""), sep=";")
# # # # # 	    
# # # # # #	    x <- match(fromVolcano$GeneID, allgeneprobes[,"gene"])
# # # # # #	    ExprData <- ExprData[which(is.na(x) == FALSE),]
# # # # # #	    Pvalues <- allPvalues[which(is.na(x) == FALSE)]
# # # # # #	    coeffs <- allcoeffs[which(is.na(x) == FALSE)]
# # # # # 	    #coeffs <- fromVolcano$coefficients[which(is.na(x) == FALSE)]
# # # # # #	    geneprobes <- allgeneprobes[x[which(!is.na(x))],]
# # # # # 	    
# # # # # #	    assign(esetName, new("ExpressionSet",exprs=as.matrix(ExprData),phenoData=new("AnnotatedDataFrame",data=pData)))
# # # # #             assign(esetName4save, new("ExpressionSet", exprs=as.matrix(ExprData), 
# # # # # 				      phenoData=new("AnnotatedDataFrame",data=pData),
# # # # # 				      featureData=new("AnnotatedDataFrame",
# # # # # # 					    data=as.data.frame(rownames(ExprData), row.names=rownames(ExprData), stringsAsFactors=FALSE))
# # # # # 					    data=fData)
# # # # # 				      ))
# # # # # 
# # # # # 	    resultdf <- try(plot_with_properties_eset(eset=get(esetName4save), geneprobes, prop, web=FALSE, dataType,
# # # # # 						filterchk, basement, includeORexclude, ifCluster, ifLR, exp, db=platform, path_to_write,
# # # # # 						dataFormat, up, down, window.size, diffGenesSelected=TRUE, Pvalues,
# # # # # 						eset4save=get(esetName4save), allPvalues=Pvalues, allcoeffs=coeffs,
# # # # # 						allgeneprobes=geneprobes, MLinfo=MLinfo))
# # # # # 	    endA <- Sys.time()
# # # # # 	    print(paste("start time for plot_with_properties_eset",startA))
# # # # # 	    print(paste("stop time for plot_with_properties_eset",endA))
# # # # # 	    print(paste("total time for plot_with_properties_eset",endA-startA))
# # # # # 	    print(paste("done",resultdf$Remark[1],"at time",Sys.time()))
# # # # # 
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- paste(length(unique(geneprobes[,"gene"])))
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "NoOfGenes"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- paste(length(geneprobes[,"gene"]))
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "NoOfProbes"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- paste(as.vector(coeffs),collapse="DELIM")
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficients"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- paste(geneprobes[match(names(coeffs),geneprobes[,"probe_name"]),"gene"],collapse="DELIM")
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficientsGeneID"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- paste(geneprobes[match(names(coeffs),geneprobes[,"probe_name"]),"symbol"],collapse="DELIM")
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "coefficientsGeneName"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- levels(status)[1] 
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "level1"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]]<- levels(status)[2] 
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "level2"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- length(status[status==resultdf[1,"level1"]]) 
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "level1no"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- length(status[status==resultdf[1,"level2"]]) 
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "level2no"
# # # # # 	    
# # # # # 	    resultdf[[(length(resultdf)+1)]] <- paste(as.vector(Pvalues),collapse="DELIM")
# # # # # 	    colnames(resultdf)[dim(resultdf)[2]] <- "Pvals"
# # # # # 	    print(paste("Current sessionID is: ", resultdf$sessionID))
# # # # # # 	    file.copy(paste(path_to_write, "genelist_temp_volcano_", previousID, ".xls", sep=""), paste0(path_to_write, "/temp_volcano_", resultdf$sessionID, ".xls"))
# # # # # ##	    file.copy(paste(path_to_write, "/temp_volcano_", previousID, ".xls", sep=""), paste0(path_to_write, "/temp_volcano_", resultdf$sessionID, ".xls"))
# # # # # # 	    write.table(fromVolcano, file=paste0(path_to_write, "/temp_volcano_", resultdf$sessionID, ".xls"), row.names=FALSE, sep="\t", quote=FALSE)
# # # # # 	    
# # # # # 	    return(resultdf)
# # # # # 	    
# # # # # 	    }
	    
	} else { # iflength(t) > 0
		print("No genes found")
                remark <- "Filtered zero"
                gpgene <- "NA"
                gpprobe <- "NA"
                dc2 <- "NA"
                pkeggvalue <- "NA"
                sample_after_filtering <- "NA"
                filteredSamples <- "NA"
                sessionID <- "NA"
                pvalueLR <- "NA"
                noOfGenes <- "NA"
                noOfProbes <- "NA"
		coeffs <- "NA"
		coeffsGeneID <- "NA"
		coeffsGeneName <- "NA"
                sample_after_filtering <- 0
                filteredSamples <- "NA"
		level1 <- "NA" 
		level2 <- "NA" 
		Pvals <- "NA"
                x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR,noOfGenes,noOfProbes,coeffs,coeffsGeneID,coeffsGeneName,level1,level2,Pvals)
                colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR","NoOfGenes","NoOfProbes","coefficients","coefficientsGeneID", "coefficientsGeneName","level1","level2","Pvals")
		print(paste("Current sessionID is: ", x$sessionID))
                return(x)
                
	}
###########################################################
###########################################################
## Starting multigroup!!!

# # # # # } else {
# # # # # 
# # # # # if(pValue=="" || is.null(pValue)) pValue <- NULL
# # # # # # if(pValue=="") pValue <- NULL
# # # # # if(foldchange=="" || is.null(foldchange)) foldchange <- NULL
# # # # # # if(foldchange=="") foldchange <- NULL
# # # # # 
# # # # # if(reAnalyse) {
# # # # #     print(paste("...reAnalysing"))
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
# # # # # #    previousID <- MLinfo$previousID
# # # # # #    reAnalyse <- MLinfo$reAnalyse
# # # # # }
# # # # # MLinfo <- list(exp=exp, prop=prop, filterchk=filterchk, basement=basement, includeORexclude=includeORexclude, ifCluster=ifCluster,
# # # # # 		up=up, down=down, window.size=window.size, path_to_write=path_to_write, pORqVal=pORqVal, 
# # # # # 		previousID=previousID, reAnalyse=reAnalyse, Multi=TRUE)
# # # # # 	## log
# # # # # 	a<-Sys.time()
# # # # # 	print(paste("master_diffGenes called with following parameters at",a))
# # # # #         print(paste("exp",exp))
# # # # #         print(paste("prop",prop))
# # # # #         print(paste("filterchk",filterchk))
# # # # #         print(paste("basement",basement))
# # # # #         print(paste("includeORexclude",includeORexclude))
# # # # #         print(paste("pValueCutoff",pValue))
# # # # #         print(paste("foldchange",foldchange))
# # # # #         print(paste("pORqVal",pORqVal))
# # # # #         print(paste("path_to_write",path_to_write))
# # # # #         print(paste("previousID",previousID))
# # # # # 	print(paste("before get pdata"))
# # # # # 
# # # # # 	startQuery <- Sys.time()
# # # # # ##	library(RMySQL)
# # # # # ##	source("http://eh3.uc.edu/r/GRS_prod.R")
# # # # # ##	source("http://eh3.uc.edu/r/mod_Rserve_prod.R")
# # # # # 	## Get dataType of the experiment
# # # # # 	#!servSet <- getServerSettings(test)
# # # # # 	#!mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# # # # # #@!	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
# # # # # 	
# # # # # #    if (!test) {
# # # # # #	mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="10.165.4.231", port = 4040, password="public")
# # # # # #    } else {
# # # # # #	mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="gimm2.ketl.uc.edu", password="public")}
# # # # #     
# # # # #         #mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="db2.ketl.uc.edu", port = 3306, password="public")
# # # # #         #!sql <- paste("select Platform,DataType,DataFormat,eset_path from ExperimentMetadata where Experiment = '",exp,"'",sep="")
# # # # # #@!        sql <- paste("select platformDb,DataType,DataFormat,eset_path from ExperimentMetadata where experimentName = '",exp,"'",sep="")
# # # # #         #!rs <- DBI::dbSendQuery(mycon, sql)
# # # # #         #!dt <- DBI::fetch(rs, n = -1)
# # # # #         dt <- expInfo(exp, test)
# # # # #         dataType <- dt["DataType"]
# # # # #         dataFormat <- dt["DataFormat"]
# # # # #         platform <- dt["Platform"]
# # # # #         esetpath <- dt["eset_path"]
# # # # #         esetpath <- try(sub("mnt", "opt", esetpath))
# # # # #         #!DBI::dbDisconnect(mycon)
# # # # # 	db <- platform
# # # # # 	## Get pData using getDataForGRS function in GRS.R
# # # # # 	pData <- getpData(exp, platform, test)
# # # # # 	rownames(pData) <- gsub(paste0(exp, "_"), "", pData$MeasurementName) ## Mehdi: Added 2018 for matching samples
# # # # # 	pData[is.na(pData)] <- "NA"
# # # # # 	
# # # # # 	print("after get pdata")
# # # # # 
# # # # # 	## Get expression data using getDataForGRS function in GRS.R
# # # # #        	esetName <- paste(exp,"_",platform,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# # # # #        	print(paste(esetName))
# # # # # 
# # # # # 	if(!is.null(glist)) {
# # # # # 		    eset <- try(getData(exp=exp, glist=glist, filterchk=filterchk, includeORexclude=includeORexclude, homoloGenes=homoloGenes, debuging=debuging))
# # # # # # 		    if(length(grep("^Error.*",eset))){  # if eset not found at the specified location
# # # # # 		    if(class(eset) != "ExpressionSet"){  # if eset not found at the specified location
# # # # # 			    esetpath=""  # query from database
# # # # # 		    } else {
# # # # # # 			    eset <- as(get(eset),"ExpressionSet")
# # # # # # 			    if(sum(is.na(match(rownames(pData),colnames(Biobase::exprs(eset))))) > 0)  # not matching
# # # # # # 			    {
# # # # # # 				    colnames(Biobase::exprs(eset)) <- rownames(pData)
# # # # # # 			    }
# # # # # # 			    Biobase::pData(eset) <- pData[colnames(Biobase::exprs(eset)), , drop=FALSE]
# # # # # # 			    fdat <- annotateFeatures(exp, debuging=debuging)
# # # # # # 			    fdat <- fdat[match(rownames(Biobase::exprs(eset)), fdat$PROBE), ]
# # # # # # 			    rownames(fdat) <- fdat$PROBE
# # # # # # 			    Biobase::fData(eset) <- fdat
# # # # # 			    assign(esetName,eset)
# # # # # 			    probes <- rownames(Biobase::exprs(eset))
# # # # # 			    geneprobes <- get_gene_probe_eset(db, unique(probes),test)
# # # # # 		    }
# # # # # 	} else {
# # # # # 	    if(esetpath != "")  # esetpath present. directly load the eset from the directory
# # # # # 	    {
# # # # # 		eset <- try(load(esetpath))
# # # # # 		if(length(grep("^Error.*",eset))){  # if eset not found at the specified location
# # # # # 			esetpath=""  # query from database
# # # # # 		} else {
# # # # # 			eset <- as(get(eset),"ExpressionSet")
# # # # # 			if(sum(is.na(match(rownames(pData),colnames(Biobase::exprs(eset))))) > 0)  # not matching
# # # # # 			{
# # # # # # 				colnames(Biobase::exprs(eset)) <- rownames(pData)
# # # # # 				Biobase::exprs(eset) <- Biobase::exprs(eset)[ , rownames(pData), drop=FALSE]
# # # # # 			}
# # # # # 			Biobase::pData(eset) <- pData[colnames(Biobase::exprs(eset)), , drop=FALSE]
# # # # # 			fdat <- annotateFeatures(exp, debuging=debuging)
# # # # # # 			fdat <- fdat[match(rownames(Biobase::exprs(eset)), fdat$PROBE), ]
# # # # # 			fdat <- fdat[rownames(Biobase::exprs(eset)), ]
# # # # # # 			rownames(fdat) <- fdat$PROBE
# # # # # 			Biobase::fData(eset) <- fdat
# # # # # 			assign(esetName,eset)
# # # # #                 	probes <- rownames(Biobase::exprs(eset))
# # # # # 			geneprobes <- get_gene_probe_eset(db, unique(probes), test)
# # # # # 		}
# # # # # 	    }
# # # # # 	    if(esetpath == "")  # eset path info not provided then obtain the data by querying the database
# # # # # 	    {
# # # # # 		eData <- getDataForGRS(exp,filterchk,basement,includeORexclude,samples=rownames(pData), platform)
# # # # # 		eData <- eData[which(eData[,"Gene"] != "" ),]
# # # # # 
# # # # # 		## Create eset
# # # # # #		ord <- order(pData[,prop])
# # # # # #       	pData <- pData[ord,]
# # # # # #        	NoOfProbes <- dim(eData)[1]/dim(pData)[1]
# # # # # #       	geneNames <- eData[1:NoOfProbes,"Gene"]
# # # # # 
# # # # # 		ExprData <- NULL
# # # # #         	ExprData <- matrix(eData[,"Value"],ncol=dim(pData)[1])
# # # # # 		# unqVals <- as.character(unique(eData[,1]))
# # # # # 		# unqIndex <- match(unqVals,rownames(pData))
# # # # # 
# # # # #         	rownames(ExprData) <- eData[1:length(unique(eData[,"Probe"])), "Probe"] 
# # # # # 		colnames(ExprData) <- rownames(pData)
# # # # #         	ExprData<-as.data.frame(ExprData)
# # # # # 
# # # # # 		## Call functions from plot_with_properties_eset wrapper from mod_Rserve.R
# # # # # 		probes <- rownames(ExprData)
# # # # # 		geneprobes <- get_gene_probe_eset(db,unique(probes), test)
# # # # # 		ExprData <- ExprData[geneprobes[,"probe_name"],]
# # # # # 
# # # # # # 		library(Biobase)
# # # # # # 		loadNamespace("Biobase")
# # # # #         	assign(esetName,new("ExpressionSet",exprs=as.matrix(ExprData),
# # # # # 				    phenoData=new("AnnotatedDataFrame",data=pData),
# # # # # 				    featureData=new("AnnotatedDataFrame",data=as.data.frame(rownames(ExprData), row.names=rownames(ExprData), stringsAsFactors=FALSE))))
# # # # # 
# # # # # 		endQuery <- Sys.time()
# # # # # 		print(paste("start time for query and eset",startQuery))
# # # # # 		print(paste("end time for query and eset",endQuery))
# # # # # 		print(paste("total time for query and eset",endQuery - startQuery))
# # # # # 	    #save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
# # # # # 	    } # if no esetpath found
# # # # # 	}
# # # # # 
# # # # # ## # # # # geneprobes4volcano <- geneprobes
# # # # # #geneprobes4volcano <- geneprobes4volcano
# # # # # 	## filter eset according to criterion and then run anova
# # # # # 	f <- filter_eset(get(esetName),filterchk,basement,includeORexclude)
# # # # #  	
# # # # # #  	remove NA  or missing data
# # # # #  	ExprData <- Biobase::exprs(f)[rownames(geneprobes), , drop=FALSE]
# # # # # 	countingMissingData <- apply(ExprData,1,function(x) sum(is.na(x)))
# # # # # 	nonMissingData <- (countingMissingData != length(colnames(ExprData)))
# # # # # 	ExprData <- ExprData[nonMissingData, ,drop=FALSE]
# # # # # 	ExprData <- as.matrix(ExprData)
# # # # # #	rownames(ExprData) <- as.vector(rownames(ExprData)[names(nonMissingData)])
# # # # # 
# # # # #        	pData <- Biobase::pData(f)
# # # # # 	pData[is.na(pData)] <- "NA"
# # # # #        	fData <- Biobase::fData(f)[rownames(ExprData), , drop=FALSE]
# # # # # 	fData <- as.data.frame(apply(fData, 2, as.character), stringsAsFactors=FALSE, row.names=rownames(fData))
# # # # # # # # # # 	    if (reAnalyse) { ## need to check
# # # # # # # # # # 		print(paste("...reAnalysing"))
# # # # # # # # # # 		#fromVolcano <- read.csv(paste(path_to_write, "genelist_multi_volcano_", previousID, ".csv", sep=""), sep=";")
# # # # # # # # # # 		#geneprobes <- geneprobes[match(fromVolcano$Probes, geneprobes$probe_name),]
# # # # # # # # # # 		#ExprData <- ExprData[geneprobes[!is.na(geneprobes$probe_name),"probe_name"], ,drop=FALSE] #Biobase::exprs(f)[geneprobes[,"probe_name"],]
# # # # # # # # # # 		#ExprData <- ExprData[unique(rownames(ExprData)), ,drop=FALSE]
# # # # # # # # # # 	    }
# # # # # 
# # # # # #  	remove NA  or missing data
# # # # # 	
# # # # # 	print("before Anova")
# # # # # 	## Step 1 compute anova
# # # # # 	startanova <- Sys.time()
# # # # # ##	library(limma)
# # # # # ##	library(genefilter)
# # # # # 	status <- as.factor(as.character(pData[,prop]))
# # # # #         #design <- model.matrix(~ -1+status)
# # # # #         design <- model.matrix(~status)
# # # # #         fit <- limma::lmFit(ExprData, design = design)
# # # # #  #      junk<-summary(lm(ExprData[1,]~status))
# # # # #  #	junk<-summary(lm(ExprData[1,]~-1+status))
# # # # #         fit2 <- NULL
# # # # # 	try(fit2 <- genefilter::rowFtests(ExprData,status))
# # # # # #       try(fit2 <- eBayes(fit),TRUE)
# # # # # 
# # # # # 	if(is.null(fit2))
# # # # #         {
# # # # #                 #return(paste("Error in anova",geterrmessage()))
# # # # # 		print("Error in anova")
# # # # # 		remark <- "Error in anova"
# # # # # 		gpgene <- "NA"
# # # # # 		gpprobe <- "NA"
# # # # # 		dc2 <- "NA"
# # # # # 		pkeggvalue <- "NA"
# # # # # 		sample_after_filtering <- "NA"
# # # # # 		filteredSamples <- "NA"
# # # # # 		sessionID <- "NA"
# # # # # 		pvalueLR <- "NA"
# # # # # 		noOfGenes <- "NA"
# # # # # 		noOfProbes <- "NA"
# # # # #                 sample_after_filtering <- 0
# # # # #                 filteredSamples <- "NA"
# # # # #                 x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR,noOfGenes,noOfProbes)
# # # # #                 colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR","NoOfGenes","NoOfProbes")
# # # # # 		print(paste("Current sessionID is: ", x$sessionID))
# # # # # 
# # # # #                 return(x)
# # # # # 
# # # # #         }
# # # # #         
# # # # # 	fit2 <- fit2[complete.cases(fit2), ,drop=FALSE]
# # # # # 	fit <- fit$coefficients[rownames(fit2), ,drop=FALSE]
# # # # # #fit <- fit
# # # # # #fit2 <- fit2
# # # # # #ExprData <- ExprData
# # # # # 
# # # # # 	print("after Anova")
# # # # #         ## if exactly 2 levels use p.values else use F.p.value
# # # # # # # # # #         if(length(unique(pData[,prop])) == 2) {
# # # # # # # # # #                 #Pvalues <- fit2$p.value[,2]
# # # # # # # # # #                 Pvalues <- fit2$p.value
# # # # # # # # # #                 names(Pvalues)<-rownames(fit2)
# # # # # # # # # #         } else {
# # # # #                 #Pvalues <- fit2$F.p.value
# # # # #                 #names(Pvalues)<-names(fit2$p.value[,2])
# # # # #                 Pvalues <- fit2$p.value
# # # # #                 names(Pvalues)<-rownames(fit2)
# # # # # # # # # #         }
# # # # # 
# # # # # 	if(pORqVal=="q") {
# # # # # ##		library(qvalue)
# # # # # 		Pvalues <- qvalue::qvalue(p=Pvalues,lambda=0)$qvalues
# # # # # 	}
# # # # # #	Pvalues <- Pvalues[!is.na(Pvalues)]
# # # # # 	
# # # # # 	## Creating full signature!
# # # # # 	Pvalues <- sort(Pvalues[!is.na(Pvalues)])
# # # # # # # # # # 	geneprobes <- geneprobes[!is.na(geneprobes[,"probe_name"]), ,drop=FALSE] ## Mehdi: it was ...ed up. NA.s always showed up as a gene
# # # # # # # # # # 	geneprobes <- geneprobes[!is.na(geneprobes[,"gene"]), ,drop=FALSE] ## 
# # # # # # # # # # 	Pvalues <- Pvalues[rownames(geneprobes)]
# # # # # 	geneprobes <- geneprobes[names(Pvalues), ]
# # # # # geneprobes4volcano <- geneprobes
# # # # # 	
# # # # #         multi_volcano <- geneprobes
# # # # #         multi_volcano <- cbind(multi_volcano, 
# # # # # 				"Pvals" = Pvalues, 
# # # # # 				"coefficients" = apply(fit[names(Pvalues),-1, drop=FALSE], 1, function(x) max(x, na.rm=TRUE))
# # # # # 				)
# # # # #         colnames(multi_volcano) <- c("ID_geneid", "PROBE", "Name_GeneSymbol", "Description", "Significance_pvalue", "Value_LogDiffExp")
# # # # #         # [c("probe_name", ), ]
# # # # #         
# # # # # #         multi_volcano <- data.frame("Probes" = rownames(fit2), #[match(rownames(ExprData), names(Pvalues))], 
# # # # # # 				    "Pvals" = fit2$p.value, #[match(rownames(ExprData),names(Pvalues))], 
# # # # # # 				    "coefficients" = apply(fit[,-1, drop=FALSE], 1, function(x) max(x, na.rm=TRUE)), stringsAsFactors=FALSE)
# # # # #         #multi_volcano <- data.frame("Probes" = names(ptvol), "Pvals" = ptvol, "coefficients" = ltvol[vol], stringsAsFactors=FALSE)
# # # # # 	#multi_volcano <- cbind(names(ptvol), ptvol, ltvol[vol])
# # # # # 	#colnames(multi_volcano) <- c("Probes", "Pvals", "coefficients")
# # # # # 	multi_volcano <- multi_volcano[!is.na(as.numeric(multi_volcano$Significance_pvalue)), ,drop=FALSE]
# # # # # 	multi_volcano <- multi_volcano[!is.na(as.numeric(multi_volcano$Value_LogDiffExp)), ,drop=FALSE]
# # # # # # # # # # 	mat <- match(multi_volcano$Probes, geneprobes4volcano$probe_name)
# # # # # # # # # # 	multi_volcano$geneID <- geneprobes4volcano[mat, "gene"]
# # # # # # # # # # 	multi_volcano$geneName <- geneprobes4volcano[mat, "symbol"]
# # # # # # # # # # 	multi_volcano$Description <- geneprobes4volcano[mat, "description"]
# # # # # # # # # # 	multi_volcano <- multi_volcano[, c("Probes", "geneID", "geneName", "coefficients", "Pvals", "Description")]
# # # # # # # # # # 	multi_volcano <- multi_volcano[match(geneprobes4volcano$probe_name, multi_volcano$Probes), ]
# # # # # # # # # # 	colnames(multi_volcano) <- c("PROBE", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue", "Description")
# # # # # #@!	multi_volcano <- multi_volcano[, c(2, 1, 3, 6, 5, 4)]
# # # # # 	
# # # # # #multi_volcano <- multi_volcano
# # # # # # # # # # 	volcano_name <- paste(path_to_write,"multi_volcano_", as.character(resultdf$sessionID), ".xls", sep = "")
# # # # # # # # # # 	write.table(multi_volcano, file=volcano_name, row.names=FALSE, sep="\t", quote=FALSE)
# # # # # 	## End creating full signature should be saved later!!!
# # # # # 
# # # # # 	
# # # # # 	## Step 2 filter accd to p/q value and/or fold change
# # # # # 	# Filter by user defined pvalue
# # # # # 	##Mehdi: changed to use only the first 100 genes
# # # # # 	if (!is.null(pValue)) {
# # # # # 	    Pvalues <- Pvalues[which(Pvalues <= as.numeric(pValue))]
# # # # # # 	    ptvol <- Pvalues
# # # # # 	    print(paste("filtering Pvals")) 
# # # # # 	} else {
# # # # # # 	    pt <- ptvol <- Pvalues
# # # # # 	    print(paste("Pvals null"))
# # # # # #		if (is.null(foldchange))  pt <- sort(Pvalues)[1: min(100, length(Pvalues))]; print(paste("first 100"))
# # # # # #		if (!is.null(foldchange))  pt <- Pvalues; print(paste("Ps null, Fs not!"))
# # # # # 	}
# # # # # 	# Filter by user defined foldchange
# # # # # 	##Mehdi: removed foldchange from multigroup pipeline
# # # # #         if (!is.null(foldchange)) {
# # # # # 	    print(paste("filtering Folds"))
# # # # # 	    if(length(unique(pData[,prop])) == 2) {
# # # # # 		    #lt<-fit2$coefficients[which(abs(fit2$coefficients[,2])>log2(as.numeric(foldchange))),]
# # # # # 		    lt <- fit[which(abs(fit[,2])>=(as.numeric(foldchange))), , drop=FALSE] ## Mehdi: I removed log2 from as.numeric(foldchange)
# # # # # 		    ltvol <- fit[, 2, drop=FALSE]
# # # # # 	
# # # # # 	    } else {  # if more than 2 levels take the mx fold change
# # # # # 		    #lt <- fit2$coefficients[which(abs(apply(fit2$coefficients[,-1],1,function(x) max(marray::rm.na(x)))) > log2(as.numeric(foldchange))),]
# # # # # #		    lt <- fit$coefficients[which(abs(apply(fit$coefficients[,-1, drop=FALSE],1,function(x) max(marray::rm.na(x)))) >= log2(as.numeric(foldchange))), , drop=FALSE]
# # # # # 		    lt <- fit[which(abs(apply(fit[,-1, drop=FALSE],1,function(x) max(x, na.rm=TRUE))) >= (as.numeric(foldchange))), , drop=FALSE] ## Mehdi: I removed log2 from as.numeric(foldchange)
# # # # # 		    ltvol <- fit# apply(fit$coefficients[,-1, drop=FALSE],1,function(x) max(x, na.rm=TRUE))
# # # # # 	    }
# # # # # 	} else {
# # # # # 	    lt <- ltvol <- fit
# # # # # 	    print(paste("Folds null"))
# # # # # 	}
# # # # # #	lt <- lt[complete.cases(lt), ,drop=FALSE]
# # # # # 	
# # # # # #ptvol <- ptvol
# # # # # #ltvo <- ltvol
# # # # # 	endanova <- Sys.time()
# # # # # 	print(paste("start time for anova",startanova))
# # # # # 	print(paste("end time for anova",endanova))
# # # # # 	print(paste("total time for anova",endanova - startanova))
# # # # # 
# # # # # 	## Intersect
# # # # # # 	if (is.null(pValue) & is.null(foldchange) & is.null(glist)) {
# # # # # # 		pt <- sort(Pvalues)[1: min(100, length(Pvalues))]
# # # # # # 		print(paste("first 100 probes"))
# # # # # # 		t <- names(pt) 
# # # # # # #	    } else if (is.null(pValue) & !is.null(foldchange)) { t <- rownames(lt) 
# # # # # # #	    } else if (!is.null(pValue) & is.null(foldchange)) { t <- names(pt) 
# # # # # # 	    } else { t <- intersect(names(pt), rownames(lt)) } # Mehdi: this was the original statement
# # # # # 
# # # # # 	## Intersect
# # # # # 	
# # # # # 	if (is.null(pValue) & is.null(foldchange) & is.null(glist)) {
# # # # # 		pt <- Pvalues[1: min(100, length(Pvalues))]
# # # # # 		print(paste("first 100 probes"))
# # # # # 		t <- names(pt) 
# # # # # 	} else {
# # # # # 		t <- intersect(names(Pvalues), rownames(lt))
# # # # # 	} # Mehdi: this was the original statement
# # # # # 	
# # # # # 	geneprobes <- geneprobes[t, ,drop=FALSE]
# # # # # 	
# # # # # 	    
# # # # # 	if(length(t) > 0) {
# # # # # #tt <- t
# # # # # 	# Get genes
# # # # # # 	geneprobes <- geneprobes[match(t,geneprobes[,"probe_name"]), ,drop=FALSE]
# # # # # # 	geneprobes <- geneprobes[!is.na(geneprobes[,"probe_name"]), ,drop=FALSE] ## Mehdi: it was ...ed up. NA.s always showed up as a gene
# # # # # # # # # # 	x <- match(geneprobes[,"probe_name"], rownames(ExprData))
# # # # # # # # # # 	if(sum(is.na(x)) > 0) {
# # # # # # # # # # 		x <- x[which(is.na(x) == FALSE)]
# # # # # # # # # # 	#	colnames(ExprData) <- paste(exp,"_",colnames(ExprData),sep="" )
# # # # # # # # # # 
# # # # # # # # # # 	}
# # # # # # # # # # 	ExprData <- ExprData[x, , drop=FALSE]
# # # # # 	
# # # # # 	ExprData <- ExprData[t, , drop=FALSE]
# # # # # 	fData <- fData[match(rownames(ExprData), rownames(fData)), ]
# # # # # 
# # # # # 	# Create eset
# # # # #         esetName <- paste(exp,"_", platform, "_", paste(unlist(strsplit(gsub(":", " ", date()), " ")), collapse="_", sep=""), sep="")
# # # # # # 	esetName4save <- paste(esetName, "4save", sep="")
# # # # # 	assign(esetName,new("ExpressionSet", exprs=as.matrix(ExprData), 
# # # # # 			    phenoData=new("AnnotatedDataFrame",data=pData),
# # # # # 			    featureData=new("AnnotatedDataFrame",
# # # # # # 			    data=as.data.frame(rownames(ExprData), row.names=rownames(ExprData), stringsAsFactors=FALSE))
# # # # # 			    data=fData)
# # # # # 			    ))
# # # # # #       assign(esetName4save, new("ExpressionSet",exprs=as.matrix(ExprData),phenoData=new("AnnotatedDataFrame",data=pData)))
# # # # #         
# # # # # 	## Use above eset in treeview, cluster,KEGG, LR etc
# # # # # 	
# # # # # 	# Do clustering only if number of genes is less than 1000
# # # # # 	if(length(unique(geneprobes[,"gene"])) > 1000) {
# # # # # 		ifCluster <- "none"}
# # # # # 	if(dim(ExprData)[1] == 1) {
# # # # # 		ifCluster <- "col"}
# # # # # 	
# # # # # 	# LR not applicable
# # # # # 	ifLR <- 0
# # # # # 	# pass the computed Pvalues
# # # # # 	Pvalues <- Pvalues[match(rownames(ExprData),names(Pvalues))]
# # # # # 	if (is.null(pValue)) { pValue <- 1 }
# # # # # 	
# # # # # #Pvalues <- Pvalues
# # # # # #geneprobes <- geneprobes
# # # # # #ExprData <- ExprData
# # # # # 	startA <- Sys.time()
# # # # # 	resultdf <- try(plot_with_properties_eset(eset=get(esetName),geneprobes,prop,web=FALSE,dataType, 
# # # # # 						filterchk,basement,includeORexclude,ifCluster,ifLR,exp,db=platform,path_to_write,
# # # # # 						dataFormat,up,down,window.size,diffGenesSelected=TRUE,Pvalues,
# # # # # 						PvalueCutoff=as.numeric(pValue), MLinfo=MLinfo, eset4save=get(esetName)))
# # # # # 	endA <- Sys.time()
# # # # # 	print(paste("start time for plot_with_properties_eset",startA))
# # # # # 	print(paste("stop time for plot_with_properties_eset",endA))
# # # # # 	print(paste("total time for plot_with_properties_eset",endA-startA))
# # # # # 	print(paste("done",resultdf$Remark[1],"at time",Sys.time()))
# # # # # #	resultdf <- cbind(resultdf, c(length(unique(geneprobes[,"gene"]))))
# # # # # 	resultdf <- cbind(resultdf, c(length(unique(geneprobes[,"gene"]))))
# # # # #         colnames(resultdf)[dim(resultdf)[2]] <- "NoOfGenes"
# # # # # 	resultdf <- cbind(resultdf, c(length(t)))
# # # # #         colnames(resultdf)[dim(resultdf)[2]] <- "NoOfProbes"
# # # # #         
# # # # # # # # # #         ## Mehdi: to create a volcano
# # # # # # # # # #         print(paste("Current sessionID is: ", resultdf$sessionID))
# # # # # # # # # #         #vol <- match(names(ptvol), names(ltvol))
# # # # # # # # # #         #vol <- ltvol[!is.na(match(names(ltvol),names(ptvol)))]
# # # # # # # # # #         multi_volcano <- data.frame("Probes" = rownames(fit2), #[match(rownames(ExprData), names(Pvalues))], 
# # # # # # # # # # 				    "Pvals" = fit2$p.value, #[match(rownames(ExprData),names(Pvalues))], 
# # # # # # # # # # 				    "coefficients" = apply(fit[,-1, drop=FALSE], 1, function(x) max(x, na.rm=TRUE)), stringsAsFactors=FALSE)
# # # # # # # # # #         #multi_volcano <- data.frame("Probes" = names(ptvol), "Pvals" = ptvol, "coefficients" = ltvol[vol], stringsAsFactors=FALSE)
# # # # # # # # # # 	#multi_volcano <- cbind(names(ptvol), ptvol, ltvol[vol])
# # # # # # # # # # 	#colnames(multi_volcano) <- c("Probes", "Pvals", "coefficients")
# # # # # # # # # # 	multi_volcano <- multi_volcano[!is.na(as.numeric(multi_volcano$Pvals)), ,drop=FALSE]
# # # # # # # # # # 	multi_volcano <- multi_volcano[!is.na(as.numeric(multi_volcano$coefficients)), ,drop=FALSE]
# # # # # # # # # # 	mat <- match(multi_volcano$Probes, geneprobes4volcano$probe_name)
# # # # # # # # # # 	multi_volcano$geneID <- geneprobes4volcano[mat, "gene"]
# # # # # # # # # # 	multi_volcano$geneName <- geneprobes4volcano[mat, "symbol"]
# # # # # # # # # # 	multi_volcano$Description <- geneprobes4volcano[mat, "description"]
# # # # # # # # # # 	multi_volcano <- multi_volcano[, c("Probes", "geneID", "geneName", "coefficients", "Pvals", "Description")]
# # # # # # # # # # 	multi_volcano <- multi_volcano[match(geneprobes4volcano$probe_name, multi_volcano$Probes), ]
# # # # # # # # # # 	colnames(multi_volcano) <- c("PROBE", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue", "Description")
# # # # # # # # # # #multi_volcano <- multi_volcano
# # # # # # # # # # 	volcano_name <- paste(path_to_write,"multi_volcano_", as.character(resultdf$sessionID), ".xls", sep = "")
# # # # # # # # # # 	write.table(multi_volcano, file=volcano_name, row.names=FALSE, sep="\t", quote=FALSE)
# # # # # # # # # # 
# # # # # 	volcano_name <- paste(path_to_write,"multi_volcano_", as.character(resultdf$sessionID), ".xls", sep = "")
# # # # # 	write.table(multi_volcano, file=volcano_name, row.names=FALSE, sep="\t", quote=FALSE)
# # # # # 
# # # # #         return(resultdf)
# # # # # 
# # # # # 	} else { # iflength(t) > 0
# # # # # 		print("No genes found")
# # # # #                 remark <- "Filtered zero"
# # # # #                 gpgene <- "NA"
# # # # #                 gpprobe <- "NA"
# # # # #                 dc2 <- "NA"
# # # # #                 pkeggvalue <- "NA"
# # # # #                 sample_after_filtering <- "NA"
# # # # #                 filteredSamples <- "NA"
# # # # #                 sessionID <- "NA"
# # # # #                 pvalueLR <- "NA"
# # # # #                 noOfGenes <- "NA"
# # # # #                 noOfProbes <- "NA"
# # # # #                 sample_after_filtering <- 0
# # # # #                 filteredSamples <- "NA"
# # # # #                 x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR,noOfGenes,noOfProbes)
# # # # #                 colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR","NoOfGenes","NoOfProbes")
# # # # # 		print(paste("Current sessionID is: ", x$sessionID))
# # # # # 
# # # # #                 return(x)
# # # # # 
# # # # # 	}
# # # # # 
# # # # # }
}


