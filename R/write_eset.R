
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
#' @examples
#' res <- batchTregGRS(....)

# write_eset <- function(genes,exp,path_to_write, debuging=FALSE)
# {
# ##	library(RMySQL)
# 
# 	## Get the database name for given experiment first
# 	servSet <- getServerSettings(test)
#  	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# #	if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
#         sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#         rs <- DBI::dbSendQuery(mycon, sql)
#         tmp <- DBI::fetch(rs,n=-1)
#         DBI::dbDisconnect(mycon)
# 	db <- tmp[1,"Platform"]
# 	dataFormat <- tmp[,"DataFormat"]
# 	if(is.null(db))
# 	{
# 		return("Experiment not found. Please check the name and try again.")
# 	}
# 
# 
#         genes <- gsub(",NA","",genes)
#         genes <- gsub("NA,","",genes)
# 	geneid <- unlist(strsplit(genes,","))
# 	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)
# 	probes <- getProbes(geneid,exp,db, test)
# 	sampledata <- get_sample_properties(db,exp,meas_limits, test)
# 	parsedSData <- parse_sample_properties(propertyTags=sampledata)
# 	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
# 	rownames(featuredata) <- rownames(parsedSData)
# 	rownames(probes) <- probes$Feature_Name
# 	probes <- probes[colnames(featuredata), ]
# 	
# # 	ID_geneid <- probes[, "Gene"]
# 	idName <- geneid2symbol(paste(probes$Gene, collapse=","))
# 	idName <- idName[match(probes$Gene, idName$GeneID), ]
# # 	PROBE <- colnames(featuredata)
# 	fdata <- data.frame(ID_geneid=idName$GeneID, Name_GeneSymbol=idName$Symbol, PROBE=probes$Feature_Name , DESCRIPTION=NA , stringsAsFactors=FALSE)
# 	rownames(fdata) <- fdata$PROBE
# ##	library(Biobase)
# # 	loadNamespace("Biobase")
# 	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# 	esetName <- gsub("-", "_", esetName)
# 	assign(esetName, new("ExpressionSet",exprs=t(featuredata),phenoData=new("AnnotatedDataFrame",data=parsedSData)))
# 	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
# 	return(esetName)
# 
# }

#############
# write_eset <- function(genes,exp,path_to_write, debuging=FALSE)
# {
# ##	library(RMySQL)
# 	## Get the database name for given experiment first
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# #	if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
#         sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#         rs <- DBI::dbSendQuery(mycon, sql)
#         tmp <- DBI::fetch(rs,n=-1)
#         DBI::dbDisconnect(mycon)
# 	db <- tmp[1,"Platform"]
# 	dataFormat <- tmp[,"DataFormat"]
# 	if(is.null(db))
# 	{
# 		return("Experiment not found. Please check the name and try again.")
# 	}
# 
#         genes <- gsub(",NA","",genes)
#         genes <- gsub("NA,","",genes)
# # 	geneid <- unlist(strsplit(genes,","))
# 	geneid <- genes
# 	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)
# 	probes <- getProbes(geneid,exp,db, test)
# 	sampledata <- get_sample_properties(db,exp,meas_limits, test)
# 	parsedSData <- parse_sample_properties(propertyTags=sampledata)
# 	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
# 	rownames(featuredata) <- rownames(parsedSData)
# 	rownames(probes) <- probes$Feature_Name
# 	probes <- probes[colnames(featuredata), ]
# 	
# # 	ID_geneid <- probes[, "Gene"]
# 	idName <- geneid2symbol(paste(probes$Gene, collapse=","))
# 	idName <- idName[match(probes$Gene, idName$GeneID), ]
# # 	PROBE <- colnames(featuredata)
# 	fdata <- data.frame(ID_geneid=idName$GeneID, Name_GeneSymbol=idName$Symbol, PROBE=probes$Feature_Name , DESCRIPTION=NA , stringsAsFactors=FALSE)
# 	rownames(fdata) <- fdata$PROBE
# ##	library(Biobase)
# 	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# 	esetName <- gsub("-", "_", esetName)
# 	assign(esetName, new("ExpressionSet", exprs=t(featuredata), phenoData=new("AnnotatedDataFrame", data=parsedSData), 
# 			      featureData=new("AnnotatedDataFrame", data=fdata)))
# 	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
# 	esetToGct1.3(get(esetName), paste0("/filteredeset_", esetName), path_to_write)
# 	return(esetName)
# 
# }
##########


########## working 12-1-2017
# write_eset <- function(genes,exp,path_to_write, debuging=FALSE)
# {
# ##	library(RMySQL)
# 	loadNamespace("Biobase")
# 	## Get the database name for given experiment first
# # 	servSet <- getServerSettings(test)
# # 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# # #	if (!test) {
# # #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# # #	else {
# # #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
# #         sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
# #         rs <- DBI::dbSendQuery(mycon, sql)
# #         tmp <- DBI::fetch(rs,n=-1)
# #         DBI::dbDisconnect(mycon)
# 	tmp <- expInfo(exp, test)
# 	db <- tmp$Platform
# 	dataFormat <- tmp$DataFormat
# 	if(is.null(db))
# 	{
# 		return("Experiment not found. Please check the name and try again.")
# 	}
# 
#         genes <- gsub(",NA","",genes)
#         genes <- gsub("NA,","",genes)
# # 	geneid <- unlist(strsplit(genes,","))
# 	geneid <- genes
# 	
#     if(dataFormat=="MaxD") {
# # 	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)
# 
#         ### getProbes function
# 	servSet <- getServerSettings(test) 
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# 	pi <- DBI::dbGetQuery(mycon, paste0("SELECT ID FROM Platforms WHERE Platform = '", db, "'"))[1,1]
# 	probes <- DBI::dbGetQuery(mycon, paste0("SELECT Probe, Gene FROM ProbeGene WHERE Platform_ID = '", pi, "' AND Gene IN (",geneid,")"))
#         colnames(probes) <- c("Feature_Name","Gene")
#         DBI::dbDisconnect(mycon)
#         ###
# 
# 
# # 	probes <- getProbes(geneid,exp,db, test)
# # 	sampledata <- get_sample_properties(db,exp,meas_limits, test)
# # 	parsedSData <- parse_sample_properties(propertyTags=sampledata)
# # 	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
# # 	rownames(featuredata) <- rownames(parsedSData)
# # 	rownames(probes) <- probes$Feature_Name
# # 	probes <- probes[colnames(featuredata), ]
# 	if(nrow(probes)==0)
# 	{
# 		return("None of the genes found.")
# 	}
# 	eset <- get(load(sub("/mnt", "/opt", tmp$eset_path))) # , envir=environment()
# 	edat <- Biobase::exprs(eset)
# 	edat <- edat[probes$Feature_Name, , drop=FALSE]
# 	pdat <- Biobase::pData(eset)[,, drop=FALSE]
# # 	fdat <- Biobase::fData(eset)
# 	
# # 	ID_geneid <- probes[, "Gene"]
# 	idName <- geneid2symbol(paste(probes$Gene, collapse=","))
# 	idName <- idName[match(probes$Gene, idName$GeneID), ]
# # 	PROBE <- colnames(featuredata)
# 	fdata <- data.frame(ID_geneid=idName$GeneID, Name_GeneSymbol=idName$Symbol, PROBE=probes$Feature_Name , DESCRIPTION=NA , stringsAsFactors=FALSE)
# 	rownames(fdata) <- fdata$PROBE
# ##	library(Biobase)
# 	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# 	esetName <- gsub("-", "_", esetName)
# # 	assign(esetName, new("ExpressionSet", exprs=t(featuredata), phenoData=new("AnnotatedDataFrame", data=parsedSData), 
# # 			      featureData=new("AnnotatedDataFrame", data=fdata)))
# 	assign(esetName, new("ExpressionSet", exprs=edat, phenoData=new("AnnotatedDataFrame", data=pdat), 
# 			      featureData=new("AnnotatedDataFrame", data=fdata)))
# 	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
# 	esetToGct1.3(get(esetName), paste0("/filteredeset_", esetName), path_to_write)
# 	return(esetName)
# 	
# 	
#     } else {
# 	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)
# 	probes <- getProbes(geneid,exp,db, test)
# 	sampledata <- get_sample_properties(db,exp,meas_limits, test)
# 	parsedSData <- parse_sample_properties(propertyTags=sampledata)
# 	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
# 	rownames(featuredata) <- rownames(parsedSData)
# 	rownames(probes) <- probes$Feature_Name
# 	probes <- probes[colnames(featuredata), ]
# 	
# # 	ID_geneid <- probes[, "Gene"]
# 	idName <- geneid2symbol(paste(probes$Gene, collapse=","))
# 	idName <- idName[match(probes$Gene, idName$GeneID), ]
# # 	PROBE <- colnames(featuredata)
# 	fdata <- data.frame(ID_geneid=idName$GeneID, Name_GeneSymbol=idName$Symbol, PROBE=probes$Feature_Name , DESCRIPTION=NA , stringsAsFactors=FALSE)
# 	rownames(fdata) <- fdata$PROBE
# ##	library(Biobase)
# 	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# 	esetName <- gsub("-", "_", esetName)
# 	assign(esetName, new("ExpressionSet", exprs=t(featuredata), phenoData=new("AnnotatedDataFrame", data=parsedSData), 
# 			      featureData=new("AnnotatedDataFrame", data=fdata)))
# 	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
# 	esetToGct1.3(get(esetName), paste0("/filteredeset_", esetName), path_to_write)
# 	return(esetName)
#     }
# 
# }
##########

##########
write_eset <- function(genes,exp,path_to_write, homoloGenes=TRUE, debuging=FALSE)
{
##	library(RMySQL)
# 	loadNamespace("Biobase")
	## Get the database name for given experiment first
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# #	if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
#         sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#         rs <- DBI::dbSendQuery(mycon, sql)
#         tmp <- DBI::fetch(rs,n=-1)
#         DBI::dbDisconnect(mycon)
	tmp <- expInfo(exp, debuging)
	db <- tmp$Platform
	dataFormat <- tmp$DataFormat
	if(is.null(db))
	{
		return("Experiment not found. Please check the name and try again.")
	}

        genes <- gsub(",NA","",genes)
        genes <- gsub("NA,","",genes)
# 	geneid <- unlist(strsplit(genes,","))
	geneid <- genes
	
    if(dataFormat=="MaxD") {
# 	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)

        ### getProbes function
	servSet <- getServerSettings(debuging) 
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#@!	pi <- DBI::dbGetQuery(mycon, paste0("SELECT ID FROM Platforms WHERE Platform = '", db, "'"))[1,1]
#@!	probes <- DBI::dbGetQuery(mycon, paste0("SELECT Probe, Gene FROM ProbeGene WHERE Platform_ID = '", pi, "' AND Gene IN (",geneid,")"))
	probes <- DBI::dbGetQuery(mycon, paste0("SELECT distinct Probe, Gene FROM ProbeGene PG,Platforms P where P.Platform ='", db,"' AND PG.Platform_ID = P.ID AND Gene in (", geneid,") AND Gene NOT LIKE \"%/%\""))
	
        colnames(probes) <- c("Feature_Name","Gene")
        probes <- probes[!is.na(as.integer(probes$Gene)), ] ## Dec,2017
        DBI::dbDisconnect(mycon)
        ###


# 	probes <- getProbes(geneid,exp,db, test)
# 	sampledata <- get_sample_properties(db,exp,meas_limits, test)
# 	parsedSData <- parse_sample_properties(propertyTags=sampledata)
# 	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
# 	rownames(featuredata) <- rownames(parsedSData)
# 	rownames(probes) <- probes$Feature_Name
# 	probes <- probes[colnames(featuredata), ]
	if(nrow(probes)==0)
	{
		return("None of the genes found.")
	}
	
# 	eset <- get(load(sub("/mnt", "/opt", tmp$eset_path))) # , envir=environment()
# 	edat <- Biobase::exprs(eset)
# 	edat <- edat[probes$Feature_Name, , drop=FALSE]
# 	pdat <- Biobase::pData(eset)[,, drop=FALSE]
# # 	fdat <- Biobase::fData(eset)
# 	
# # 	ID_geneid <- probes[, "Gene"]
# 	idName <- geneid2symbol(paste(probes$Gene, collapse=","), description=TRUE)
# 	idName <- idName[match(probes$Gene, idName$GeneID), ]
# # 	PROBE <- colnames(featuredata)
# 	fdata <- data.frame(ID_geneid=idName$GeneID, Name_GeneSymbol=idName$Symbol, PROBE=probes$Feature_Name , DESCRIPTION=idName$Description , stringsAsFactors=FALSE)
# 	rownames(fdata) <- fdata$PROBE
# ##	library(Biobase)
# 	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# 	esetName <- gsub("-", "_", esetName)
# # 	assign(esetName, new("ExpressionSet", exprs=t(featuredata), phenoData=new("AnnotatedDataFrame", data=parsedSData), 
# # 			      featureData=new("AnnotatedDataFrame", data=fdata)))
# 	assign(esetName, new("ExpressionSet", exprs=edat, phenoData=new("AnnotatedDataFrame", data=pdat), 
# 			      featureData=new("AnnotatedDataFrame", data=fdata)))
	
	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
	esetName <- gsub("-", "_", esetName)
	assign(esetName, getData(exp=exp, glist=genes, homoloGenes=homoloGenes, debuging=debuging))
	
	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
	esetToGct1.3(get(esetName), paste0("/filteredeset_", esetName), path_to_write)
	return(esetName)
	
	
    } else {
	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)
	probes <- getProbes(geneid,exp,db, debuging)
	sampledata <- get_sample_properties(db,exp,meas_limits, debuging)
	parsedSData <- parse_sample_properties(propertyTags=sampledata)
	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
	rownames(featuredata) <- rownames(parsedSData)
	rownames(probes) <- probes$Feature_Name
	probes <- probes[colnames(featuredata), ]
	
# 	ID_geneid <- probes[, "Gene"]
	idName <- geneid2symbol(paste(probes$Gene, collapse=","), debuging=debuging)
	idName <- idName[match(probes$Gene, idName$GeneID), ]
# 	PROBE <- colnames(featuredata)
	fdata <- data.frame(ID_geneid=idName$GeneID, Name_GeneSymbol=idName$Symbol, PROBE=probes$Feature_Name , DESCRIPTION=NA , stringsAsFactors=FALSE)
	rownames(fdata) <- fdata$PROBE
##	library(Biobase)
	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
	esetName <- gsub("-", "_", esetName)
	assign(esetName, new("ExpressionSet", exprs=t(featuredata), phenoData=new("AnnotatedDataFrame", data=parsedSData), 
			      featureData=new("AnnotatedDataFrame", data=fdata)))
	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
	esetToGct1.3(get(esetName), paste0("/filteredeset_", esetName), path_to_write)
	return(esetName)
    }

}

# write_eset <- function(genes,exp,path_to_write, debuging=FALSE)
# {
# ##	library(RMySQL)
# 
# 	## Get the database name for given experiment first
# 	servSet <- getServerSettings(test)
# #@! 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
# #	if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
# 
# 	tmp <- DBI::dbGetQuery(mycon, paste0("select platformDb,DataFormat from ExperimentMetadata where experimentName = '",exp,"'"))
#         DBI::dbDisconnect(mycon)
# 	db <- tmp[1,"platformDb"]
# 	dataFormat <- tmp[,"DataFormat"]
# 	if(is.null(db))
# 	{
# 		return("Experiment not found. Please check the name and try again.")
# 	}
# 
#         genes <- gsub(",NA","",genes)
#         genes <- gsub("NA,","",genes)
# 	geneid <- unlist(strsplit(genes,","))
# 	meas_limits <- initiate(db,exp,web=FALSE,dataFormat, debuging=debuging)
# 	probes <- getProbes(geneid,exp,db, test)
# 	sampledata <- get_sample_properties(db,exp,meas_limits, test)
# 	parsedSData <- parse_sample_properties(propertyTags=sampledata)
# 	featuredata <- feature_data(db,exp,feature=probes[,1:2],meas=parsedSData[,"MeasurementIDs"],meas_limits=NULL, debuging=debuging)
# 	rownames(featuredata) <- rownames(parsedSData)
# ##	library(Biobase)
# 	esetName <- paste(exp,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
# 	esetName <- gsub("-", "_", esetName)
# 	assign(esetName, new("ExpressionSet",exprs=t(featuredata),phenoData=new("AnnotatedDataFrame",data=parsedSData)))
# 	save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
# 	return(esetName);
# 
# }
