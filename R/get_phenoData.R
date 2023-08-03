
#@' A function for ...
#@'
#@' This function allows you to ...
#@' @param exp This is basically the experiment name.
#@' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#@'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#@' @param authors M. Fazel-Najafabadi
#@' @keywords filter
#@' @export
#@' @examples
#@' res <- get_phenoData(....)

#@get_phenoData <- function(exp, debuging=FALSE) {
  
#@servset <- getServerSettings(test)
#@  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = servset$port, host=servset$host, password='public')#@
#@  sql <- paste("select eset_path from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#@  eset_path <- DBI::dbGetQuery(mycon, sql)[1,"eset_path"]

#@  esetpath <- gsub("/mnt/","/opt/",eset_path)
#@  esetname <- load(esetpath)
#@  eset <- get(esetname)
#@  dat <- Biobase::pData(eset)
#@  dat[is.na(dat)] <- "NA"
#@  res <- list(header = colnames(dat), rows = dat)  #works single header
#@  return(res)
#@}

#' A function for ...
#'
#' This function allows you to ...
#' @param genelist A character string containing a list of Entrez gene IDs separated by commas.
#' @param experiment The name of selected experiment from portal as a single character string.
#' @param path_to_write This parameter specifies the path which user wants to save the results in.
#' @param up,down Specifies the maximum number of up and down regulated genes.
#' @param window.size ...
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

# get_phenoData <- function(genelist="none", experiment, path_to_write="/mnt/raid/tmp", up=4000, down=1000, window.size=50, debuging=FALSE)
# {
# 	print(paste("Preprocess called with following parameters at",Sys.time()))
# 	print(paste("genelist",genelist))
# 	print(paste("experiment",experiment))
# 	print(paste("path_to_write",path_to_write))
# 	print(paste("up",up))
# 	print(paste("down",down))
# 	print(paste("window.size",window.size))
# 
# ##	library(RMySQL)
# 	requireNamespace("Biobase", quietly = TRUE)
# 	## First get homologenes
# 	genes <- find_genes_hgenes_in_platform(exp=experiment, glist=genelist, debuging=debuging)
# 
# 	if(is.null(genes))
# 	{
# 		return(list(c("No genes found","NA","NA"),NULL))
# 	}
# 	#print(genes)
#         
# 	## Get the database name for given experiment first
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# #        if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
#         sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",experiment,"'",sep="")
#         rs <- DBI::dbSendQuery(mycon, sql)
#         tmp <- DBI::fetch(rs,n=-1)
#         DBI::dbDisconnect(mycon)
# 
#         db <- tmp[1,"Platform"]
# 	dataFormat <- tmp[1,"DataFormat"]
#         if(is.null(db))
#         {
# 		return(list(c("Experiment not found. Please check the name and try again.","NA","NA"),NULL))
#         }
# 	
# 	if(!is.null(genes))
# 	{
# 		if(dataFormat == "MaxD")
# 		{
# 			esetName <- write_eset(genes,experiment,path_to_write, test)
# 		} else {
# 			esetName <- create_gff_eset(glist=genes,experiment,up,down,window.size,path_to_write, test)
# 		}
# 	}
# 	if(esetName == "Experiment not found. Please check the name and try again.")
# 	{
# 		return(list(c("Experiment not found. Please check the name and try again.","NA","NA"),NULL))
# 
# 	}
# 
# 	if(esetName == "No genes found")
# 	{
# 		return(list(c("No genes found","NA","NA"),NULL))
# 
# 	}
# 	val <- c(esetName,db,genes)
# 	#print(val)
# 	## Get the property values for given experiment
#         propVals <- NULL
# 	load(paste(path_to_write, esetName, ".RData", sep=""))
# 	eset <- get(esetName)
# 	for(p in 1:dim(Biobase::pData(eset))[2])
# 	{
# 		 propVals <- rbind(propVals,c(colnames(Biobase::pData(eset))[p],paste(unique(Biobase::pData(eset)[,p]),collapse=",",sep="") ))
# 	}
# 	colnames(propVals) <- c("Property","Unique values")
# 	returnList <- list(val,propVals)
# 	return(returnList)
# 	
# }
# 






# get_phenoData_old <- function(exp, genelist="none", path_to_write="/mnt/raid/tmp/", up=4000, down=1000, window.size=50, debuging=FALSE)
# {
# 	experiment <- exp
# 	print(paste("get_phenoData called with following parameters at",Sys.time()))
# 	print(paste("genelist",genelist))
# 	print(paste("experiment",experiment))
# 	print(paste("path_to_write",path_to_write))
# 	print(paste("up",up))
# 	print(paste("down",down))
# 	print(paste("window.size",window.size))
# 
# ##	library(RMySQL)
# 	requireNamespace("Biobase", quietly = TRUE)
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# 	sql <- paste("select Organism,Platform,DataFormat,eset_path from ExperimentMetadata where Experiment = '",experiment,"'",sep="")
# 	rs <- DBI::dbSendQuery(mycon, sql)
# 	tmp <- DBI::fetch(rs,n=-1)
# # 	tmp <- tmp
# #	DBI::dbDisconnect(mycon)
# 	dataFormat <- tmp[1,"DataFormat"]
# 	org <- tmp[1,"Organism"]
# 	platform <- tmp[1,"Platform"]
# 	eset_path <- tmp[1, "eset_path"]
# 	eset_path <- gsub("/mnt/","/opt/",eset_path)
# 	print(paste("dataFormat", dataFormat))
# 	print(paste("Organism", org))
# 	print(paste("Platform", platform))
# 	print(paste("eset_path", eset_path))
# 
# 	## First get homologenes
# 	if(genelist != "none") { # !is.null(genelist) & 
# 	    genes <- find_genes_hgenes_in_platform(exp=experiment, glist=genelist, debuging=debuging)
# 	} else { #if (eset_path=="") {#if(genelist=="none") {#genes <- "none"
# 	t1 <- Sys.time()
# 	tmp2 <- NULL
# 	sql2 <- paste("SELECT distinct Gene FROM ProbeGene PG,Platforms P where P.Platform ='",platform,"' AND PG.Platform_ID = P.ID AND Gene NOT LIKE '%/%'",sep="")
# 	tmp2 <- DBI::dbGetQuery(mycon, sql2)
# #@	platformgenes <- paste(tmp2[,1],collapse=",",sep="")
# #@	genes <- gsub(",,",",",platformgenes)
# #@	genes <- gsub(",,",",",platformgenes)
# 	platformgenes <- tmp2[,1]
# 	genes <- na.omit(platformgenes)
# 	genes <- sample(genes, min(20, length(tmp2[,1])))
# 	print(paste("getting genes in:", Sys.time()-t1))
# # 	genes <- "none"
# 	print(paste("no genes", length(unique(tmp2[,1]))))
# # 	print(paste(genes))
# 	} # else genes <- "none"
# #  	tmp2 <- tmp2
# 
# 	if(is.null(genes))
# 	{
# 		return(list(c("No genes found","NA","NA"),NULL))
# 	}
# 	#print(genes)
#         
# 	## Get the database name for given experiment first
# #@	servSet <- getServerSettings(test)
# #@	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# #        if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')}
# #	else {
# #	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", host="gimm2.ketl.uc.edu", password='public')}
# #@        sql <- paste("select Platform,DataFormat from ExperimentMetadata where Experiment = '",experiment,"'",sep="")
# #@        rs <- DBI::dbSendQuery(mycon, sql)
# #@        tmp <- DBI::fetch(rs,n=-1)
#         DBI::dbDisconnect(mycon)
# 
#         db <- platform
# #@	dataFormat <- tmp[1,"DataFormat"]
#         if(is.null(db))
#         {
# 		return(list(c("Experiment not found. Please check the name and try again.","NA","NA"),NULL))
#         }
# 	
# 	if(!is.null(genes))# & genes != "none")
# 	{
# 		if(dataFormat == "MaxD")
# 		{
# 		    if(genelist == "none" & eset_path != "")
# 		    {
# 			esetName <- sub(".RData", "", unlist(strsplit(eset_path, "/"))[length(unlist(strsplit(eset_path, "/")))])
# 		    } else {
# 			esetName <- write_eset(genes,experiment,path_to_write, test)
# 		    }
# 		} else {
# 			esetName <- create_gff_eset(glist=genes,experiment,up,down,window.size,path_to_write, test)
# 		}
# 	}
# 	
# 	print(paste("esetName", esetName))
# 	if(esetName == "Experiment not found. Please check the name and try again.")
# 	{
# 		return(list(c("Experiment not found. Please check the name and try again.","NA","NA"),NULL))
# 
# 	}
# 
# 	if(esetName == "No genes found")
# 	{
# 		return(list(c("No genes found","NA","NA"),NULL))
# 
# 	}
# #@	val <- c(esetName,db,genes)
# 	#print(val)
# 	## Get the property values for given experiment
# #@        propVals <- NULL
# 	if (genelist == "none" & eset_path != ""){ # dataFormat == "MaxD") {
# 	    tmp <- load(eset_path)
# 	    eset <- get(tmp)
# 	} else { 
# 	    tmp <- load(paste0(path_to_write, esetName, ".RData")) }
# 	    eset <- get(tmp)
# #@	for(p in 1:dim(Biobase::pData(eset))[2])
# #@	{
# #@		 propVals <- rbind(propVals,c(colnames(Biobase::pData(eset))[p],paste(unique(Biobase::pData(eset)[,p]),collapse=",",sep="") ))
# #@	}
# 	
# 	dat <- Biobase::pData(eset)
# 	dat[is.na(dat)] <- "NA"
# 	res <- list(header = colnames(dat), rows = unique(dat))  #works single header
# 	return(res)
# 
# #@	colnames(propVals) <- c("Property","Unique values")
# #@	returnList <- list(val,propVals)
# #@	return(returnList)
# 	
# }



get_phenoData <- function(exp, debuging=FALSE) {

## 	servSet <- getServerSettings(test)
## 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
# # # 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
##	sql <- paste("select PropertyName, PropertyValue from MeasurementMetadata where Experiment = '",exp,"'",sep="")## #@!
# # #@!	sql <- paste("select Measurement_ID, PropertyName, PropertyValue from MeasurementMetadata where Experiment = '",exp,"'",sep="")
##	rs <- DBI::dbGetQuery(mycon, sql)# # #@!
# # 	dataformat <- DBI::dbGetQuery(mycon, paste0("select DataFormat FROM ExperimentMetadata WHERE Experiment = '",exp,"'"))[1,1]
# # 	if (dataformat == "MaxD") {
# # 	    rs <- DBI::dbGetQuery(mycon, paste0("select Measurement_ID, PropertyName, PropertyValue FROM MeasurementMetadata WHERE Experiment = '",exp,"'"))
# # #	tmp <- tmp
##	tmp <- split(rs, rs$PropertyName)## #@!
##	for (i in 1:length(tmp)) colnames(tmp[[i]])[2] <- names(tmp)[i]## #@!
##	tmp <- data.frame(sapply(tmp, function(j) j[,2]), stringsAsFactors=F)## #@!
	tmp <- getpData(exp, debuging)
# # 	
# # 	    rs <- as.data.frame(lapply(reshape::cast(rs, Measurement_ID ~ PropertyName, value="PropertyValue")[-1], as.vector), stringsAsFactors=FALSE)
# # 	  } else {
# # 	      requireNamespace("Biobase", quietly = TRUE)
# # 	      db <- DBI::dbGetQuery(mycon, paste0("select Platform FROM ExperimentMetadata WHERE Experiment = '",exp,"'"))[1,1]
# # 	      mycon2 <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
# # # 	      somegenes <- DBI::dbGetQuery(mycon2, "SELECT F.Name, F.ID, G.Name FROM Gene G, GeneList GL, Feature F, Reporter R WHERE G.Name NOT LIKE \"%/%\" AND GL.Gene_ID = G.ID AND R.Gene_List_ID = GL.ID AND F.Reporter_ID = R.ID")
# # 	      somegenes <- DBI::dbGetQuery(mycon2, "SELECT GeneID FROM refGene")
# # 	      somegenes <- paste(sample(na.omit(as.integer(somegenes[,1])), 100), collapse=",")
# # 	      tmp <- create_gff_eset(somegenes,exp,up=4000,down=1000,window.size=50,path_to_write="/mnt/raid/tmp",debuging=debuging, write=FALSE)
# # 	      rs <- Biobase::pData(tmp)
# # 	      rs <- unique(as.data.frame(rs, stringsAsFactors=F))
# # 	      row.names(rs) <- NULL
# # 	  }
	tmp[is.na(tmp)] <- "<NA>"
 	res <- list(header=colnames(tmp), rows=tmp)## #@!
# # 	res <- list(header=colnames(rs), rows=rs)
}

#################
# getProbes <- function(geneid=NULL,experiment=NULL,db=NULL, debuging=FALSE){
#         probes <- NULL
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
# #        if (!test) {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="10.165.4.231", port = 4040, password='public')}
# #	else {
# #	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')}
# 	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
# ## Mehdi: This function was a bit funny! ;-)
# #         sql <- paste0("SELECT PG.Probe, PG.Feature_ID, PG.Gene FROM ProbeGene PG, ExperimentMetadata EM, PlatformInfo PI WHERE EM.experimentName = '", experiment, "' AND PI.Platform = EM.PlatformDb AND PG.Platform_ID = PI.ID AND Gene IN (",geneid,") AND Gene NOT LIKE \"%/%\"")
# #         sql <- paste0("SELECT PG.Probe, PG.Feature_ID, PG.Gene FROM ProbeGene PG, ExperimentMetadata EM, PlatformInfo PI WHERE Gene IN (",geneid,") AND Gene NOT LIKE \"%/%\" AND PI.Platform = EM.PlatformDb AND PG.Platform_ID = PI.ID AND EM.experimentName = '", experiment, "'")
# # Mehdi: using these queries takes about 2.5 seconds, no good!
# 
# 	pn <- DBI::dbGetQuery(mycon, paste0("SELECT platformDb FROM ExperimentMetadata WHERE experimentName = '", experiment, "'"))[1,1]
# 	pi <- DBI::dbGetQuery(mycon, paste0("SELECT ID FROM PlatformInfo WHERE Platform = '", pn, "'"))[1,1]
# 	probes <- DBI::dbGetQuery(mycon, paste0("SELECT Probe, Feature_ID, Gene FROM ProbeGene WHERE Platform_ID = '", pi, "' AND Gene IN (",geneid,")"))
#         colnames(probes) <- c("Feature_Name","Feature_ID","Gene")
#         DBI::dbDisconnect(mycon)
#         return(probes)
# }
# 
# 
