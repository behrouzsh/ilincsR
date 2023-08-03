
#' A function for ...
#'
#' This function allows you to ...
#' @param db db is the database name in the GenomicsPortals.
#' @param exp This is basically the experiment name.
#' @param prop is the property selected from the experiment.
#' @param filterchk filterchk is in the form of a string. It can have multiple property of the pData to filter
#'	based on. Each pair of property:value should be separated by a comma ",". Pairs are saparated
#'	by colon ":". An example ot filterchk is :
#'	filterchk="property1:value1,property1:value2,property2:value1"
#' @param basement basement Levels of comparison, not relevant here!
#' @param ifCluster This parameter defines in which way the data should be clustered.
#' @param includeORexclude This argument basically is designed to filter the pData based on what 
#'	filterchk is or selected property of eset when is set to "1", or the reverse selection of filterchk 
#'	when set to "2". 
#'	It also can have NULL or "n" values when there is no property selected, in this case 
#'	the function will return the original ExpressionSet.
#' @param path_to_write This parameter specifies the path which user wants to save the results.
#' @param up,down Specifies the maximum number of up and down regulated genes.
#' @param window.size ...
#' @param display Which resalts to show. It should be a vector of column indices or "all" for everything.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param authors M. Fazel-Najafabadi
#' @keywords filter
#' @export
#' @examples
#' res <- analyze(exp="TCGA_KIRC_RNASeqV2", filterchk="n", prop="n", ifCluster="both", includeORexclude="2", 
#'			path_to_write= "/mnt/raid/tmp/", esetName="TCGA_BRCA_RNASeqV2_tcgaIlluminaHiSeq_Wed_Sep_13_12_30_57_2017", display=c(1,6,10))


analyze <- function(db=NULL, exp, prop="n",filterchk=NULL, basement=NULL, ifCluster="both", includeORexclude="1",
		    esetName, ifLR="0", path_to_write, up=4000, down=1000, window.size=50, display="all", debuging=FALSE)
{
	if(!is.null(filterchk)) filterchk <- gsub("<NA>", "NA", filterchk)
	a <- Sys.time()
	print(paste("Analyze called with following parameters at",a))
# 	print(paste("db",db))
	print(paste("exp",exp))
	print(paste("path_to_write",path_to_write))
	print(paste("up",up))
	print(paste("down",down))
	print(paste("window.size",window.size))
	print(paste("prop",prop))
	print(paste("filterchk",filterchk))
	print(paste("ifCluster",ifCluster))
	print(paste("includeORexclude",includeORexclude))
	print(paste("esetName",esetName))
	print(paste("ifLR",ifLR))

##	library(RMySQL)
	tmpeset <- load(paste(path_to_write,esetName,".RData",sep=""))
# 	servSet <- getServerSettings(test)
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#         sql <- paste("select Platform,DataType,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#         rs <- DBI::dbSendQuery(mycon, sql)
#         dt <- DBI::fetch(rs, n = -1)
#         dataType <- dt[1,"DataType"]
# 	dataFormat <- dt[1,"DataFormat"]
# 	db <- dt[1,"Platform"]
#         DBI::dbDisconnect(mycon)
        tmp <- expInfo(exp, debuging)
        dataType <- tmp$DataType
        dataFormat <- tmp$DataFormat
        db <- tmp$Platform
        
	esetString <- esetName
	esetName <- get(tmpeset)
#@!	meas_limits <- initiate(db,exp,web=TRUE,dataFormat, debuging=debuging)

	totalSamples <- NULL
	filteredSamples <- NULL
	if (dataFormat == "Gff") ifCluster <- "row"
	totalSamples <- dim(Biobase::exprs(esetName))[2]
	probes <- rownames(Biobase::exprs(esetName))
	geneprobes <- get_gene_probe_eset(db, probes, debuging)
	
#@!	geneprobes <- Biobase::fData(esetName)
#@!	colnames(geneprobes) <- c("gene","symbol","probe_name","description")
	resultdf <- plot_with_properties_eset(esetName, geneprobes, prop,web=FALSE, dataType, filterchk, basement, includeORexclude,
					    ifCluster, ifLR, exp, db, path_to_write, dataFormat, up, down, window.size)
	
	if(resultdf[1,"Remark"] != "Filterd zero" & resultdf[1,"Remark"] != "all missing")
	{
		samples <- resultdf[1,"FilteredSamples"]
		newEset <- write_fitered_eset(esetString, samples, db, exp, path_to_write)	
		resultdf <- cbind(resultdf, newEset)
		f <- filter_eset(eset=esetName, filterchk=filterchk, basement=basement, includeORexclude=includeORexclude, 
				 new.esetName=paste0("filteredeset_", resultdf$sessionID), path_to_write="/mnt/raid/tmp/")
		means <- apply(data.matrix(Biobase::exprs(f)),1,median,na.rm=T)
		centeredData <- sweep(data.matrix(Biobase::exprs(f)),1,means,"-")
		Biobase::exprs(f) <- centeredData

		save(f, file=paste(path_to_write, "/filteredeset_", resultdf$sessionID,".RData", sep=""))
		esetToGct1.3(f, paste0("/filteredeset_", resultdf$sessionID), path_to_write)
	}
	
	
# 	resultdf <- master_for_servlet(db,exp,prop,filterchk,basement,ifCluster,includeORexclude,esetName,ifLR,path_to_write,up,down,window.size, test)
	b <- Sys.time()
	print(paste("analyze done in",b-a))
	if(display != "all") return(resultdf[, display, drop=FALSE])
	return(resultdf)
	
}

# master_for_servlet <- function(db,exp,prop,filterchk,basement,ifCluster,includeORexclude,esetName, 
# 				ifLR,path_to_write,up,down,window.size, debuging=FALSE)
# {
# ##	source("http://eh3.uc.edu/r/gimmHeat.R")
# 	
# 	## Get dataType of the experiment
# 
# 		
# 	## If property is specified
# 	## Mehdi: one line for Gff clustering
# 	return(resultdf)
# }

