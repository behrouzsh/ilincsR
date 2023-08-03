
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
#' @keywords test
#' @export 
#' @examples
#' res <- batchTregGRS(....)

create_gff_eset <- function(glist,exp,up,down,window.size, path_to_write, debuging=FALSE, write=TRUE, smpls="all") {

	genes <- glist
	experiment <- exp
	db <- NULL
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname='QueryTables', host='10.165.4.231', port = 4040, password='public')}
#	else {
#	  DBI::dbConnect(RMySQL::MySQL(), user='public', dbname='QueryTables', host='gimm2.ketl.uc.edu', password='public')}
	
	sql <- paste("select Platform from ExperimentMetadata where Experiment='",experiment,"'",sep="")
# 	rs <- DBI::dbSendQuery(mycon, sql)
# 	db <- DBI::fetch(rs,n=-1)[1,1]
	db <- DBI::dbGetQuery(mycon, sql)[1,1]
	DBI::dbDisconnect(mycon)
	
	if(is.null(db)){
		print("Database not found. Check experiment name.")
		return()
	}
	
	r <-  query_gff(experiment, genelist=genes, up=up, down=down, smpls=smpls)
	if(is.null(r))
	{
		print("No genes found")
		return("No genes found")
	}
	r <- as.data.frame(r)
	
	
	cols <- "RefID,GeneID,chrom,Strand,txStart,txEnd,TSS,CHROMOSOME,chromStart,chromEnd,dataValue,startDiff,endDiff,Sample"
	colnames(r) <- unlist(strsplit(cols,","))
	us <- as.vector(unique(r[,"Sample"]))
	finalexprs <- NULL
	for(u in us)
	{
		ex <- r[which(r[,"Sample"]==u),]
		tmp <- create_expr_for_gff_sample(result=ex,up,down,window.size,missingWindow=0)
		colnames(tmp)
		finalexprs <- rbind(finalexprs,tmp)
		tmp <- NULL
	
	}
	#finalexprs
	
	
	## Get ALL the samples
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host=servSet$host, port=servSet$port, password='public')
#	mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host='10.165.4.231', port = 4040, password='public')
	sql <- paste("select S.Sample_Name from Dataset D, Sample S,Dataset_Sample DS where D.Dataset_Name = '",experiment,"' and D.Dataset_ID = DS.Dataset_ID and DS.Sample_ID = S.Sample_ID",sep="")
# 	rs <- DBI::dbSendQuery(mycon, sql)
# 	AllSamples <- DBI::fetch(rs,n=-1)
	AllSamples <- DBI::dbGetQuery(mycon, sql)
	DBI::dbDisconnect(mycon)
	if (smpls != "all") {
	  AllSamples <- AllSamples[AllSamples[,1] %in% unlist(strsplit(smpls, ",")), , drop=FALSE] 
	}
	
	cols <- seq(-up,down,by=window.size)
	cols <- cols[-length(cols)]
	colnames(finalexprs) <- c(cols,"Sample")
	cnamesFull <- NULL
	for(u in AllSamples[,1])
	{
		cnames <- paste(u,cols,sep=".")
		cnamesFull <- c(cnamesFull,cnames)
	
	}
	
	
	## Initialize matrix
	norows <- length(unique(rownames(finalexprs)))
	nocols <- (dim(finalexprs)[2]-1)* dim(AllSamples)[1]
	m <- matrix(rep(0,norows*nocols), norows, nocols)
	colnames(m) <- cnamesFull
	rownames(m) <- unique(rownames(finalexprs))
	
	## Algo to create exprs of eset
	# 2. for all probes in finalexprs
		#2.1  get sample name
		#2.2  Obtain the index of column in matrix m
		#2.3  match probe in rownames of m.
		#2.5.  insert values in m at that position
	
	
	for(t in 1:dim(finalexprs)[1])
	{
		prb <- finalexprs[t,]
		sam <- prb[length(prb)]
		colIndex <- c(which(colnames(m)==paste(sam,cols[1],sep=".")): (which(colnames(m)==paste(sam,cols[1],sep=".") ) + length(cols)-1)  )
		indexOfm <- match(rownames(finalexprs)[t],rownames(m))
		m[indexOfm, colIndex] <- as.numeric(prb[1:length(prb)-1])
	
	
	}
	
	## Create pData

	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
#	if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host='10.165.4.231', port = 4040, password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host='gimm2.ketl.uc.edu', port = 4040, password='public')}
	sql <- paste("select S.Sample_Name,SP.Property,SP.Value from Sample S, Sample_Property SP where SP.Sample_ID = S.Sample_ID",sep="")
# 	rs <- DBI::dbSendQuery(mycon, sql)
# 	props <- DBI::fetch(rs,n=-1)
	props <- DBI::dbGetQuery(mycon, sql)
	DBI::dbDisconnect(mycon)
# 	pDataCols <- unique(props[,"Property"])
	unqSamples <- unique((props[,"Sample_Name"]))
	if (smpls != "all") {
	  unqSamples <- unqSamples[unqSamples %in% unlist(strsplit(smpls, ","))]
	  props <- props[props[,"Sample_Name"] %in% unqSamples, ,drop = FALSE]
	}
	
	pdframe <- NULL
	for(u in unqSamples)
	{
		tmp <- props[props[,"Sample_Name"]==u,]
		tmpframe <- NULL
		propName <- NULL
		for(t in 1:dim(tmp)[1])
		{
			tmpframe <- cbind(tmpframe,tmp[t,"Value"])
			propName <- c(propName,tmp[t,"Property"])
			#colnames(tmpframe)<-c(colnames(tmpframe),tmp[t,"Property"])
		}
		colnames(tmpframe) <- propName
		pdframe <- rbind(pdframe,tmpframe)
	
	}
	rownames(pdframe)<-unqSamples
	
	finalpData<-NULL
	
	
	s <- colnames(m)
	t <- unlist(lapply(strsplit(s,"\\."),function(x) x[1]))
	mtch <- match(t,pdframe[,1])
	finalpData <- pdframe[t,]
	rownames(finalpData) <- s
	
	#e <- new("ExpressionSet",exprs=m,phenoData=new("AnnotatedDataFrame",data=as.data.frame(finalpData)))

	esetName <- paste(experiment,"_",db,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
	esetName <- gsub("-", "_", esetName)
	assign(esetName,new("ExpressionSet",exprs=m,phenoData=new("AnnotatedDataFrame",data=as.data.frame(finalpData))))
    if (write) {
		save(list=esetName,file=paste(path_to_write,esetName,".RData",sep=""))
		#b <- Sys.time()
		#print(b-a)
		return(esetName)
	} else {
		return(get(esetName))
	}
}
