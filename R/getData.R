
#' A function for pulling up an eset from database.
#'
#' This function allows you to recreate whole or subset of a dataset in an eset format for an experiment from database.
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
#' res <- getData(....)


getData <- function(exp, glist=NULL, filterchk=NULL, includeORexclude="1", homoloGenes=TRUE, debuging=FALSE)
{
##	library(RMySQL)

	## Get the table name for given experiment first
	if(!is.null(filterchk)) if (filterchk=="n") filterchk <- NULL
	servSet <- getServerSettings(debuging)
	esetinfo <- expInfo(exp, debuging)
 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
	tbl <- esetinfo$tblName
	platform <- esetinfo$Platform
	organism <- esetinfo$Organism
	platformID <- DBI::dbGetQuery(mycon, paste0("select ID from Platforms where Platform = '", platform, "'"))[1,1]
	mmd <- DBI::dbGetQuery(mycon, paste0("select Measurement_ID, MeasurementName, PropertyName, PropertyValue from MeasurementMetadata where Experiment = '", exp, "'"))
	if(!is.null(filterchk)) {
	    prop <- as.data.frame(t(sapply(unlist(strsplit(filterchk, ",,,")), function(i) unlist(strsplit(i, ":")))), row.names=NULL, stringsAsFactors=FALSE)
	    tp <- mmd[((mmd$PropertyName %in% prop[,1]) & (mmd$PropertyValue %in% prop[,2])), ]
	    if(as.character(includeORexclude)=="1") {
		mmd <- mmd[mmd$MeasurementName %in% tp$MeasurementName, ]
	    } else {
		mmd <- mmd[!(mmd$MeasurementName %in% tp$MeasurementName), ]
		tp <- mmd[!((mmd$PropertyName %in% prop[,1]) & (mmd$PropertyValue %in% prop[,2])), ] ## On purpose for later!
	    }
	}
	pdat <- reshape2::dcast(mmd, MeasurementName ~ PropertyName, value.var="PropertyValue")
        
	if(is.null(glist)) {
	    pg <- DBI::dbGetQuery(mycon, paste0("select Probe, Gene, Feature_ID from ProbeGene where Platform_ID = '", platformID, "'"))
	    pg <- pg[!is.na(as.integer(pg$Gene)),]
	} else {
	    glist <- gsub(",NA","", glist)
	    glist <- gsub("NA,","", glist)
	    pg <- DBI::dbGetQuery(mycon, paste0("select Probe, Gene, Feature_ID from ProbeGene where Platform_ID = '", platformID, "'", "and Gene in (", glist, ")"))
	    pg <- pg[!is.na(as.integer(pg$Gene)),]
##	    pg <- pg[!is.na(as.integer(pg$Gene)), ]
	    probeids <- paste0(pg$Feature_ID, collapse = ",")
	}
        DBI::dbDisconnect(mycon)
        if(dim(pg)[1]==0) return("Non of the genes found")
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")

	if(is.null(glist)) {
	    if(is.null(filterchk)) {
		md <- DBI::dbGetQuery(mycon, paste0("select * from ", tbl))
	    } else {
		md <- DBI::dbGetQuery(mycon, paste0("select * from ", tbl, " where Measurement_ID in (", paste(tp$Measurement_ID, collapse=","), ")"))
	    }
	    edat0 <- reshape2::dcast(md, Measurement_ID ~ Feature_ID, value.var="Value")
	    edat <- t(edat0[,-1])
	    colnames(edat) <- as.character(edat0[,1])
	} else {
	    if(is.null(filterchk)) {
		md <- DBI::dbGetQuery(mycon, paste0("select * from ", tbl, " where Feature_ID in (", probeids, ")"))
	    } else {
		md <- DBI::dbGetQuery(mycon, paste0("select * from ", tbl, " where Feature_ID in (", probeids, ") and Measurement_ID in (", paste(tp$Measurement_ID, collapse=","), ")"))
	    }
	    edat0 <- reshape2::dcast(md, Measurement_ID ~ Feature_ID, value.var="Value")
	    edat <- t(edat0[,-1])
	    rownames(edat) <- colnames(edat0)[-1]
	    colnames(edat) <- as.character(edat0[,1])
	}
	DBI::dbDisconnect(mycon)
	loadNamespace("Biobase")
	mmdd <- unique(mmd[,1:2])
	smmdd <- match(colnames(edat), as.character(mmdd$Measurement_ID))
	colnames(edat) <- mmdd$MeasurementName[smmdd]
	rownames(pdat) <- as.character(pdat$MeasurementName)
	pdat <- pdat[,-1]
	pdat <- pdat[colnames(edat), ]
	m <- match(rownames(edat), as.character(pg$Feature_ID))
	pg <- pg[m, ]
	rownames(edat) <- pg$Probe
	pg <- pg[!is.na(pg$Gene),]
	geneID <- (pg$Gene)
	if(homoloGenes) geneID <- getHomologousGenes(geneID)[, 3]
	tax <- switch(organism, human = "Hs", mouse = "Mm", rat = "Rn")
	idName <- geneid2symbol(paste(geneID, collapse=","), species=tax, debuging=debuging, description=TRUE)
	mm <- match(geneID, idName$GeneID)
	fdat <- data.frame(PROBE=pg$Probe, Name_GeneSymbol=idName$Symbol[mm], ID_geneid=as.character(geneID), DESCRIPTION=idName$Description[mm], stringsAsFactors=FALSE)
	rownames(fdat) <- fdat$PROBE
	fdat <- fdat[!is.na(fdat$ID_geneid), ] ##
	edat <- edat[rownames(fdat), , drop=FALSE]
	assign("esetName", new("ExpressionSet", 
					    exprs=as.matrix(edat), 
					    phenoData=new("AnnotatedDataFrame", data=pdat), 
					    featureData=new("AnnotatedDataFrame", data=fdat))
	)
	return(esetName)

}

