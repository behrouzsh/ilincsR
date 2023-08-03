
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
#' @examples ...
#' res <- plot_with_properties_eset(....)

plot_with_properties_eset <- function(eset, geneprobes=NULL, prop, web=FALSE, dataType, filterchk, basement, includeORexclude,
				    ifCluster, ifLR, exp, db, path_to_write, dataFormat, up, down, window.size,
				    diffGenesSelected=FALSE, Pvalues=NULL, PvalueCutoff=0.05,
				    eset4save=NULL, allPvalues=NULL, allcoeffs=NULL, allgeneprobes=NULL, MLinfo=NULL, sigDB="ilincs_sigs", chost="ilincs.org", debuging=FALSE)
{	
if (is.null(MLinfo)) MLinfo <- list(Multi=FALSE)
#@info <- list(eset=eset, geneprobes=geneprobes, prop=prop, web=web, dataType=dataType, filterchk=filterchk, basement=basement, 
#	    includeORexclude=includeORexclude, ifCluster=ifCluster, ifLR=ifLR, exp=exp, db=db, path_to_write=path_to_write, 
#	    dataFormat=dataFormat, up=up, down=down, window.size=window.size, diffGenesSelected=diffGenesSelected, 
#	    Pvalues=Pvalues, PvalueCutoff=PvalueCutoff, eset4save=eset4save, allPvalues=allPvalues, allcoeffs=allcoeffs, 
#	    allgeneprobes=allgeneprobes, MLinfo=MLinfo, debuging=debuging)
#@save(eset, geneprobes, prop, web, dataType, filterchk, basement, includeORexclude,ifCluster, ifLR, exp, db, 
#	    path_to_write, dataFormat, up, down, window.size, diffGenesSelected, Pvalues, PvalueCutoff,eset4save, 
#	    allPvalues, allcoeffs, allgeneprobes, MLinfo, test, file=paste0(path_to_write,"info.RData"))
##	library(Biobase)
#	source("http://eh3.uc.edu/r/gimmHeat.R")
	pal <- marray::maPalette(low="blue",high="yellow",mid="black")

	# Statistical analysis only if genes less than 1000
	statAnalysis <- FALSE
	if(length(unique(geneprobes[,"gene"])) <= 1000) statAnalysis <- TRUE
	print(paste("number of genes",length(unique(geneprobes[,"gene"]))))
	#####Variables to pass to kegg_pathways
	gpgene <- "NA"
	gpprobe <- "NA"
	dc2 <- "NULL"
	pkeggvalue <- "NULL"
	isPropNone <- FALSE


	#####Generating sessionID
	## Mehdi: I changed "-" with "_" and removed one "_"
        sessionID <- generateSID()

	####Filter data
	
#	expData <- Biobase::exprs(eset)
	
	f <- filter_eset(eset, filterchk, basement, includeORexclude)
# # # # # 	if (!is.null(eset4save)) {saved.eset <- filter_eset(eset4save, filterchk, basement, includeORexclude)}
	#save(f, file=paste(path_to_write, "/filteredesetinfo_", sessionID,".RData", sep="")) ## Mehdi: saving filtered eset
	filtered_eset <- Biobase::exprs(f)
	propData <- Biobase::pData(f)

	if(prop == "n") {
		pr <- apply(propData, 2, function(j) length(unique(j)))
# 		pr[pr==1] <- NA
		prp <- colnames(propData)
		prp <- prp[!pr==1]
		pr <- pr[!pr==1]
		
		
# 		prop <- colnames(propData)[1]  #assign dummy value to prop 
		prop <- prp[pr==min(pr)][1]
		isPropNone = TRUE
	}
	if(dim(filtered_eset)[2] <= 1) {
		remark <- "Filtered zero"
##Mario 		sessionID <- paste(c(strsplit(date(),split=" ")[[1]],as.integer(runif(1)*10e6)),sep="",collapse="_")
		sample_after_filtering <- 0
		filteredSamples <- "NA"
		pvalueLR <- "NA"
		x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR)
		colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR")

		return(x)

	} else {
		sample_after_filtering <- dim(filtered_eset)[2]
		filteredSamples <- paste(colnames(filtered_eset),collapse=",",sep="")
# # 		sample_in_groups <- as.data.frame(table(propData[,prop]))
# 		sample_in_groups <- paste(apply(as.data.frame(table(propData[,prop])),1,paste, collapse=",,,"), collapse="DELIM")
	}
	######### Create plotdata
	plotdata <- data.frame(colnames(filtered_eset), as.vector(propData[colnames(filtered_eset),prop]), t(filtered_eset), stringsAsFactors=FALSE)
        plotdata[,1] <- as.character(plotdata[,1])
        plotdata[,2] <- as.character(plotdata[,2])
	plotdata[is.na(plotdata[,2]),2] <- "NA"
	colnames(plotdata) <- c("SampleID",prop,rownames(filtered_eset))
#         plotdata[,2] <- as.character(plotdata[,2])


#         allColors <- c("blue","red","green","yellow","purple","brown","black","blue4","orange","green4","yellow4","aquamarine","azure4","bisque","bisque4","gold2","lightgreen","lightblue","red4","orange4","gold","gold4","firebrick","gray","gray45","khaki","pink","pink3","slategray","violet","skyblue","tomato","wheat","lightcoral","tomato3","blue","red","green","yellow","purple","brown","black","blue4","orange","green4","yellow4","aquamarine","azure4","bisque","bisque4","gold2","lightgreen","lightblue","red4","orange4","gold","gold4","firebrick","gray","gray45","khaki","pink","pink3","slategray","violet","skyblue","tomato","wheat","lightcoral","tomato3","blue","red","green","yellow","purple","brown","black","blue4","orange","green4","yellow4","aquamarine","azure4","bisque","bisque4","gold2","lightgreen","lightblue","red4","orange4","gold","gold4","firebrick","gray","gray45","khaki","pink","pink3","slategray","violet","skyblue","tomato","wheat","lightcoral","tomato3")
	data(allColors)
	nColors <- length(allColors)
	
	property_colnam <- tolower(colnames(propData))
	property_sort <- paste(prop,"_sort",sep='')
	prop_distance <- NULL
	if((dataType == "Motif Location" | dataType == "ChIP")& is.element("distance",property_colnam)) #only experiments with "distance: information
	{
		id_dist <- which(property_colnam == "distance")	
		repeat_distance <- length(grep(propData[1,id_dist],propData[,id_dist]))
		if(repeat_distance > 1) #multiple samples
		{
			prop_distance <- as.vector(propData[,id_dist])
		} else { #onlyone sample
			prop_distance <- as.vector(propData[,id_dist])
		}

		factorLevels <- unique(plotdata[,2])  # do not sort factors.
		nFactorLevels <- length(factorLevels)
	} else if(dataFormat == "Gff") {
		#from rownames of pData, get sample.distance info.
		prop_distance <- as.numeric(unlist(lapply(rownames(propData),function(x) unlist(strsplit(x,"\\."))[2])))
		factorLevels <- unique(plotdata[,2])  # do not sort factors.
		nFactorLevels <- length(factorLevels)
			
	} else if(isPropNone) { #if property is none then do not sort factors. 
		factorLevels <- unique(plotdata[,2])  # do not sort factors.
		nFactorLevels <- length(factorLevels)
	} else if(is.element(property_sort,property_colnam)) { # if the new sort order is specified
		newOrder <- as.vector(propData[colnames(filtered_eset),property_sort])
		plotdata <- plotdata[newOrder,]
		factorLevels <- unique(plotdata[,2])
		nFactorLevels <- length(factorLevels)
	} else  { #original code
		#factorLevels <- sort(unique(plotdata[,2]))
		factorLevels <- unique(plotdata[,2])
		nFactorLevels <- length(factorLevels)
		sampleOrder <- order(plotdata[,2])
		plotdata <- plotdata[sampleOrder,]
	}


	nObs <- length(plotdata[,2])
        sortedValues <- NULL
        sortedFactors <- NULL
        for(i in 1:nFactorLevels) {
                sortedValues <- c(sortedValues,plotdata[plotdata[,2]==factorLevels[i],3])
                sortedFactors <- c(sortedFactors,plotdata[plotdata[,2]==factorLevels[i],2])
        }
	nonMissing <- (!is.na(plotdata[,3]) & plotdata[,2]!="NA")

	#### factor levels check
	ifOneLevel <- FALSE
	if(length(unique(plotdata[nonMissing,2])) == 1){
		ifOneLevel <- TRUE
	}
	
	
        layoutMatrix <- matrix(c(2,1),ncol=,byrow=T)
        nObs <- dim(plotdata)[1]
        nCols <- dim(plotdata)[2]

	remark <- "Done"

	rsd <- rep("blue",nObs)
	for(i in 1:nFactorLevels){
		rsd[plotdata[,2]==factorLevels[i]] <- allColors[((i-1)%%nColors)+1]
	}
	dataClust <- data.frame(data.matrix(t(plotdata[,3:nCols])))
	countingMissingData <- apply(dataClust,1,function(x) sum(is.na(x)))
	dataClust <- cbind(colnames(plotdata)[3:nCols],colnames(plotdata)[3:nCols],dataClust)
	colnames(dataClust) <- c("ID","NAME",plotdata[,2])
	colnam <- colnames(dataClust)
	## Missing data check
	if(dataType=="ChIP" | dataType == "Motif Location")
	{
		nonMissingData <- (countingMissingData != length(colnam))
		dataClust <- dataClust[nonMissingData,]
	} else {
		conditionCheck <- (max(countingMissingData) > (nObs/2)) 
		if(conditionCheck)
		{
			nonMissingData <- (countingMissingData<((nObs)/2))
			dataClust <- dataClust[nonMissingData,]
			Pvalues <- Pvalues[nonMissingData]
			remark <- "50 percent missing"
		}
	}
	colnames(dataClust) <- colnam

	if(dim(dataClust)[1]==0)
	{
		remark <- "all missing"
		sample_after_filtering <- 0
		filteredSamples <- "NA"
		x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples)
                colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples")

		return(x)
	}
	
# # # # # 	probeInfo <- match(as.character(dataClust[,"ID"]),geneprobes[,"probe_name"])
# # # # # 	geneNames <- apply(geneprobes[na.omit(probeInfo),],1,function(x) paste(x[2],x[3],x[1],x[4],sep=" : "))
# # # # # 	dataClust <- dataClust[!is.na(probeInfo), ]
# # # # # 	dataClust[,1] <- geneprobes[na.omit(probeInfo),1]
# # # # # 	dataClust[,2] <- geneNames 
# 	probeInfo <- match(as.character(dataClust[,"ID"]),geneprobes[,"probe_name"])
	geneprobes <- geneprobes[rownames(dataClust),]
	geneNames <- apply(geneprobes,1,function(x) paste(x[2],x[3],x[1],x[4],sep=" : "))
# 	dataClust <- dataClust[!is.na(probeInfo), ]
	dataClust[,1] <- geneprobes[,1]
	dataClust[,2] <- geneNames 
	
	pvalueLR <- 5.0
	#Pvalues <- NULL
	retDf <- NULL
	## Statistical analysis..only if factor levels>1 and number of genes less than 800
	
	#if(dataType!="ChIP" & dataType != "Motif Location" & !diffGenesSelected)
	if(dataType!="ChIP" & dataType != "Motif Location" & statAnalysis)
	{ ### start for MaxD
	#	if(!ifOneLevel)
	#	{	
			retDf <- perPropertyAvg(dataClust,nFactorLevels,factorLevels,geneNames,pal,allColors,nColors,sessionID,dataType,path_to_write,Pvalues,PvalueCutoff,prop, debuging=debuging, chost=chost)## Mehdi added "prop"
			if(is.character(retDf)){
				remark <- "error in anova"
				pvalueLR <- "NA"
				ifLR <- "0"	
				Pvalues <- NULL

			} else {
				if(is.null(Pvalues)) Pvalues <- retDf[,"Pvalues"]
			}
			if(ifLR == "1" && !ifOneLevel)
			{
				
				random_eset_name <- random_genes_eset(paste(dataClust[,1],collapse=",",sep=""),exp,db,up,down,window.size,path_to_write,debuging=debuging)
#@!				dat<-list(dc=dataClust, exp=exp, db=db, path_to_write=path_to_write)
				if(random_eset_name!="No genes found")
				{

					load(paste(path_to_write,random_eset_name,".RData",sep=""))
					####Filter data
					f_random <- filter_eset(get(random_eset_name),filterchk,basement,includeORexclude)
					dataClust_random <- Biobase::exprs(f_random)
					colnames(dataClust_random) <- colnames(dataClust)[-(1:2)]
					fac <- as.factor(as.vector(colnames(dataClust_random)))
					try(Pvalues_random <- apply(dataClust_random,1,function(x) anova(lm(x~as.factor(colnames(dataClust_random))))[[5]][1]))
					#compute LR
					sampleLogPValues <- try(-log10(c(sample(Pvalues_random, length(Pvalues), replace=ifelse(length(Pvalues_random) >= length(Pvalues), FALSE, TRUE)), Pvalues))) ## Mehdi fixed for low gene number
					#sampleLogPValues <- c(sample(Pvalues_random,length(Pvalues)),Pvalues)
					sampleCat <- c(rep(0,length(Pvalues)),rep(1,length(Pvalues)))
					summary <- try(summary(glm(sampleCat~sampleLogPValues,family=binomial(link="logit"),na.action="na.omit"))) ## Mehdi added try
# 					if(length(grep(summary, "Error"))==1) {
# 						remark <- "Error in LR"
# 						x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples)
# 						colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples")
# 					    return(x)
# 					}
##					if(dim(summary$coefficients)[1] == 1)
					if(length(grep("Error", summary)) == 1)
					{
						pvalueLR <- 5.0
					} else {
						coeff <- summary$coefficients[2,1]
						pvalueLR <- summary$coefficients[2,4]

						if(coeff < 0)
						{
							pvalueLR <- 1 - pvalueLR/2
						} else {
							pvalueLR <- pvalueLR/2
						}
						pvalueLR <- signif(pvalueLR,2)
					}
				} else {
					pvalueLR <- 5.0 # NA
				
				}
			}
			
	#	} else {	
	#		remark <- "Done with one level"
	#	}

	} ### end for MaxD
	#computeLR for datatype ChIP and motive location only if sampletype is 1.

	#if(dataType=="ChIP" | dataType == "Motif Location" & !diffGenesSelected)
	if((dataType=="ChIP" | dataType == "Motif Location") & statAnalysis)
	{
		#computeLR for only if there is one level of property type.(e.g. if property type if TF, then compute LRpath for each TF separately. Instead of Pvalue take maximum peak value.
		#obtain random probes and compute maximum peak values.
		# compute LRpath
		if(nFactorLevels == 1 & ifLR == "1")
		{
			
			peakvalue <- as.vector(apply(as.matrix(dataClust[-(1:2)]),1,function(x) max(x,na.rm=TRUE)))
			random_eset_name <- random_genes_eset(paste(dataClust[,1],collapse=",",sep=""),exp,db,up,down,window.size,path_to_write)
			if(random_eset_name!="No genes found")
			{
				load(paste(path_to_write,random_eset_name,".RData",sep=""))
				####Filter data
				f_random <- filter_eset(get(random_eset_name),filterchk,basement,includeORexclude)
				dataClust_random <- Biobase::exprs(f_random)
				colnames(dataClust_random) <- colnames(dataClust)[-(1:2)]
				fac <- as.factor(as.vector(colnames(dataClust_random)))
				peakvalueRandom <- as.vector(apply(as.matrix(dataClust_random),1,function(x) max(x,na.rm=TRUE)))
				#compute LR
				sampleLogPValues <- c(peakvalue,sample(peakvalueRandom,length(peakvalue),replace=TRUE))
				sampleCat <- c(rep(1,length(peakvalue)),rep(0,length(peakvalue)))
				summary <- try(summary(glm(sampleCat~sampleLogPValues,family=binomial(link="logit")))) ## Mehdi added try
## 				if(dim(summary$coefficients)[1] == 1)
				if(length(grep("Error", summary)) == 1)
				{
					pvalueLR <- 5.0
				} else {
					coeff <- summary$coefficients[2,1]
					pvalueLR <- summary$coefficients[2,4]
		
					if(coeff < 0)
					{
						pvalueLR <- 1 - pvalueLR/2
					} else {
						pvalueLR <- pvalueLR/2
					}
					pvalueLR <- round(pvalueLR,2)
				}
			} else {
				pvalueLR <- NA
			}
		}
#save(nFactorLevels,ifLR,peakvalue,dataClust_random,peakvalueRandom,sampleLogPValues,sampleCat,summary,file="/data/srv/www/htdocs/tmp/checkLRchip.RData")
	}

	## Static heatmaps
	colnames(dataClust)[-(1:2)] <- paste(colnames(dataClust)[-(1:2)],unlist(plotdata[,1]),sep=".")
	## Mehdi: for reanalysing with a gene list
# # # # # 	if (!is.null(eset4save) || MLinfo$Multi) {
# # # # # #	previousID <- sessionID
# # # # # 	MLinfo$previousID <- sessionID
# # # # # #	save(saved.eset, previousID, dataClust, nFactorLevels, factorLevels, geneNames, pal, rsd, allColors,
# # # # # #		prop, nColors, ifCluster, prop_distance, isPropNone, path_to_write,dataFormat, allPvalues, allcoeffs, 
# # # # # #		allgeneprobes, geneprobes, file=paste(path_to_write, "/filteredesetinfo_", sessionID,".RData", sep=""))
# # # # # 
# # # # # # 	fdat <- as.data.frame(t(apply(as.matrix(dataClust[,2]), 1, function(j) unlist(strsplit(j, " : ")))), stringsAsFactors=F)
# # # # # # 	colnames(fdat) <- c("PROBE", "Name_GeneSymbol", "ID_geneid", "DESCRIPTION")
# # # # # # 	fdat <- fdat[, c(3,2,1,4)]
# # # # # # # 	fdat$ID_geneid <- as.character(fdat$ID_geneid)
# # # # # # 	rownames(fdat) <- fdat$PROBE
# # # # # # 	Biobase::fData(saved.eset) <- fdat
# # # # # 	
# # # # # # 	ed <- as.matrix(dataClust[, -(1:2)])
# # # # # # 	pd <- Biobase::pData(saved.eset)
# # # # # # 	
# # # # # # 	fd <- Biobase::fData(saved.eset)[rownames(Biobase::exprs(saved.eset)), ]
# # # # # # 	assign("saved.eset", new("ExpressionSet", exprs=,
# # # # # # 					phenoData=new("AnnotatedDataFrame",data=),
# # # # # # 					featureData=new("AnnotatedDataFrame",data=)))
# # # # # 
# # # # # ## 	means <- apply(data.matrix(dataClust[,-(1:2)]),1,median,na.rm=T)
# # # # # ## 	centeredData <- sweep(data.matrix(dataClust[,-(1:2)]),1,means,"-")
# # # # # 
# # # # # ##	Biobase::exprs(saved.eset) <- as.matrix(centeredData)
# # # # # ##	Biobase::fData(saved.eset) <- Biobase::fData(saved.eset)[rownames(Biobase::exprs(saved.eset)), ]
# # # # # ##	rownames(Biobase::pData(saved.eset)) <- paste(Biobase::pData(saved.eset)[ , prop], rownames(Biobase::pData(saved.eset)), sep=".")
# # # # # ##	Biobase::pData(saved.eset) <- Biobase::pData(saved.eset)[colnames(Biobase::exprs(saved.eset)), ]
# # # # # 	
# # # # # 	means <- apply(data.matrix(Biobase::exprs(saved.eset)),1,median,na.rm=T)
# # # # # # 	centeredData <- sweep(data.matrix(Biobase::exprs(saved.eset)),1,means,"-")
# # # # # 
# # # # # 	Biobase::exprs(saved.eset) <- as.matrix(sweep(data.matrix(Biobase::exprs(saved.eset)),1,means,"-"))
# # # # # 	Biobase::fData(saved.eset) <- Biobase::fData(saved.eset)[rownames(Biobase::exprs(saved.eset)), ]
# # # # # 	Biobase::pData(saved.eset) <- Biobase::pData(saved.eset)[colnames(Biobase::exprs(saved.eset)), ]
# # # # # 	rownames(Biobase::pData(saved.eset)) <- paste(Biobase::pData(saved.eset)[ , prop], rownames(Biobase::pData(saved.eset)), sep=".")
# # # # # 	
# # # # # 	save(saved.eset, file=paste(path_to_write, "/filteredeset_", sessionID,".RData", sep=""))
# # # # # 	esetToGct1.3(saved.eset, paste0("/filteredeset_", sessionID), path_to_write)
# # # # # 	save(MLinfo, file=paste(path_to_write, "/MLinfo_filteredeset_", sessionID,".RData", sep=""))
# # # # # 	}
	static_heatmaps(sessionID,dataClust,nFactorLevels,factorLevels,geneNames,pal,rsd,allColors,prop,nColors,ifCluster,prop_distance,isPropNone,path_to_write,dataFormat,
			sigDB=sigDB, chost=chost, debuging=debuging)

	## Download data
	writedataClust <- dataClust
	if(!(is.null(retDf)))
	{
		writedataClust <- cbind(writedataClust,retDf)
	} else if(dataFormat!="Gff" & !is.null(Pvalues)) {
		writedataClust <- cbind(writedataClust,Pvalues)
	}
	## Mehdi: changed the xls file output
	parsed <- as.data.frame(t(sapply(sapply(writedataClust$NAME, strsplit, " : "), function(i) unlist(i))))
	parsed <- parsed[, c(3, 1, 2, 4)]
# 	names(parsed) <- c("ID", "PROBE", "SYMBOL", "DESCRIPTION")
	names(parsed) <- c("ID_geneid", "PROBE", "Name_GeneSymbol", "DESCRIPTION")
#	writedataClust <- cbind(writedataClust[,1], parsed[,1], parsed[,2], parsed[,4], writedataClust[, ncol(writedataClust)-2], 
#			  writedataClust[, ncol(writedataClust)-1], writedataClust[, ncol(writedataClust)])
	xx <- ncol(writedataClust)-nFactorLevels
	writedataClust <- writedataClust[, xx:ncol(writedataClust)]
	if (nFactorLevels == 2) {writedataClust <- cbind(parsed, writedataClust[,c(1,2)], 
				  "Value_LogDiffExp"=(writedataClust[,2]-writedataClust[,1]), "Significance_pvalue"=writedataClust[, 3])
	} else {
	    writedataClust <- cbind(parsed, writedataClust)
	    colnames(writedataClust)[ncol(writedataClust)] <- "Significance_pvalue"
	    }
	write.table(writedataClust, file=paste(path_to_write, "download_data", sessionID, ".xls", sep=""), sep="\t", quote=TRUE, row.names=F, col.names=T, append=T)
# # # # # 	write.table(writedataClust, file=paste(path_to_write, "subsetSig_", sessionID, ".xls", sep=""), sep="\t", quote=FALSE, row.names=F, col.names=T, append=T)

	
	### Prepare parameters for kegg pathways
	gpgene <- paste(as.vector(dataClust[,1]),collapse="DELIM")
	gpprobe <- paste(as.vector(rownames(dataClust)),collapse="DELIM")
	dc2 <- paste(as.vector(dataClust[,2]),collapse="DELIM")		
	
	
	pkeggvalue <- rep(5,length(dataClust[,2]))

	if(dataType!="ChIP" & dataType != "Motif Location")
	{
	#	if(ifOneLevel)
	#	{			
	#		pkeggvalue <- rep(5,length(dataClust[,2]))
	#	} else {	
			pkeggvalue <- paste(as.vector(Pvalues),collapse="DELIM")
	#	}
	}

	## Mehdi: I added for volcano creation
#@	volcanoTable <- data.frame("GeneID"=dataClust[,1], "GeneNames"=geneprobes[probeInfo,3][1:dim(dataClust)[1]], "coefficients"=allcoeffs[1:length(Pvalues)], "Pvals"=Pvalues, stringsAsFactors = F)
	#gpnames <- paste(as.vector(geneprobes[probeInfo,3][1:dim(dataClust)[1]]),collapse="DELIM")
	
	#x <- data.frame(sessionID,gpgene,gpnames,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR)
	#colnames(x) <- c("sessionID","gpgene","gpnames","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR")
	x <- data.frame(sessionID,gpgene,gpprobe,dc2,pkeggvalue,remark,sample_after_filtering,filteredSamples,pvalueLR)
	colnames(x) <- c("sessionID","gpgene","gpprobe","dc2","pkeggvalue","Remark","sample_after_filtering","FilteredSamples","pvalueLR")
	
	return(x)	
		

}
