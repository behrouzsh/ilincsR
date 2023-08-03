
#' A function for "Generalized Random Set" analysis.
#'
#' This function allows you to perform the concordance analysis between two genome-scale differential expression profiles 
#'	(Query and Reference profiles). The method produces the statistical significance of the overall cocordance without 
#'	requiring specification of significance cutoffs as well as the list of genes contributing to the overall concordance. 
#'	The Query and Reference differential expression profiles can be constructed by the real-time analysis of the datasets 
#'	in Genomics Portals or can be separately uploaded in the form of the gene list with p-values of differential expression. 
#' @param data1 .
#' @param data2 .
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
#'	## not run
#'	src1 <- src2 <- "path/to/your/files" ## src2 can be set to any of available experiments on iLincs portal
#'	specie1 <- specie2 <- "human" ## or they could be different
#'	prop1 <- prop2 <- "cell_id" ## this was just a test from previous analysis, they should be grab from appropriate experiment
#'	filterchk2 <- "cell_id:MCF7,cell_id:PC3" 
#'	grs <- master_GRS("sample1.txt","sample2.txt",src1,src2,prop1,prop2,specie1,specie2,includeORexclude1=1,includeORexclude2=1,
#'			statusFile="testGRS")
#'	## end not run


master_GRS <- function(data1,data2,src1,src2,prop1,prop2,specie1,specie2,filters1=NULL,filters2=NULL,filterchk1="n",filterchk2="n", 
			includeORexclude1,includeORexclude2,statusFile, debuging=FALSE) {
##	library(CLEAN)
#	source("http://www.eh3.uc.edu/r/mod_Rserve_prod.R")
	f<-Sys.time()
	cat("1",file=paste("/data/srv/www/htdocs/tmp/",statusFile,sep=""))
	print(paste("master_GRS called with following parameters at",Sys.time()))
	print(paste("data1 ",data1))
	print(paste("data2 ",data2))
	print(paste("src1 ",src1))
	print(paste("src2 ",src2))
	print(paste("prop1 ",prop1))
	print(paste("prop2 ",prop2))
	print(paste("specie1 ",specie1))
	print(paste("specie2 ",specie2))
	print(paste("filterchk1 ",filterchk1))
	print(paste("filterchk2 ",filterchk2))
	print(paste("filters1 ",filters1))
	print(paste("filters2 ",filters2))
	print(paste("includeORexclude1 ",includeORexclude1))
	print(paste("includeORexclude2 ",includeORexclude2))
	d1Values<-NULL
	d2Values<-NULL
	d1pData<-NULL
	d2pData<-NULL
	platform1<-NULL
	platform2<-NULL
	oneTailFirst<-FALSE	
	oneTailSecond<-FALSE	
	xup<-NULL;yup<-NULL;xdn<-NULL;ydn<-NULL;yupIndex<-NULL
	levels1<-NULL;levels2<-NULL;
	## First check for one level condition
	if(src1 == "DB") 
	{
##		library(RMySQL)
		servSet <- getServerSettings(debuging)
		mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#		if (!test) {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="10.165.4.231", port = 4040, password="public")}
#		else {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="gimm2.ketl.uc.edu", password="public")}
	#	mycon <- dbConnect(MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
		sql <- paste("select Platform,eset_path from ExperimentMetadata where Experiment = '",data1,"'",sep="")
		rs <- DBI::dbSendQuery(mycon, sql)
		#platform1 <- DBI::fetch(rs,n=-1)[1,1]
		dt <- DBI::fetch(rs, n = -1)
		platform1 <- dt[1,"Platform"]
		esetpath1 <- dt[1,"eset_path"]
		DBI::dbDisconnect(mycon)

		if(is.null(platform1)) 
		{
			print("Experiment not found. Check experiment name.")
			return()
		}
		
		d1pData <- getpData(exp=data1,db=platform1)
	#	d1pData <- filter(filterchk1,filters1,parsedSData=d1pData,includeORexclude1)
		d1pData <- filter(filters1,filterchk1,parsedSData=d1pData,includeORexclude1)
		levels1 <- paste(unique(as.vector(d1pData[,prop1])),collapse=",",sep="")
		## return if one level
		if(length(unique(as.vector(d1pData[,prop1]) )) <=1)
		{
			print("Only one level")
			return("Error: Only one level for first dataset")
		}

		## check if exactly two levels
		if(length(unique(as.vector(d1pData[,prop1]) )) ==2)
		{
			print("one tailing first")
			oneTailFirst<-TRUE
		}
	}
	if(src2 == "DB") 
	{
##		library(RMySQL)
		servSet <- getServerSettings(debuging)
		mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
#		if (!test) {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="10.165.4.231", port = 4040, password="public")}
#		else {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname="QueryTables", host="gimm2.ketl.uc.edu", password="public")}
		#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")
		sql <- paste("select Platform,eset_path from ExperimentMetadata where Experiment = '",data2,"'",sep="")
		rs <- DBI::dbSendQuery(mycon, sql)
		#platform2 <- DBI::fetch(rs,n=-1)[1,1]
		dt <- DBI::fetch(rs, n = -1)
		platform2 <- dt[1,"Platform"]
		esetpath2 <- dt[1,"eset_path"]
		DBI::dbDisconnect(mycon)

		if(is.null(platform2)) 
		{
			print("Experiment not found. Check experiment name.")
			return()
		}
		
		d2pData <- getpData(exp=data2,db=platform2)
	#	d2pData <- filter(filterchk2,filters2,parsedSData=d2pData,includeORexclude2)
		d2pData <- filter(filters2,filterchk2,parsedSData=d2pData,includeORexclude2)
		levels2 <- paste(unique(as.vector(d2pData[,prop2])),collapse=",",sep="")		
		## return if one level
		if(length(unique(as.vector(d2pData[,prop2]) )) <=1)
		{
			print("Only one level")
			return("Error: Only one level for second dataset")
		}
		
		## check if exactly two levels
		if(length(unique(as.vector(d2pData[,prop2]) )) ==2)
		{
			print("one tailing second")
			oneTailSecond<-TRUE
		}
	}

	## Get data
	if(src1 == "DB") 
	{
		d1Values = NULL
		ExprData1 = NULL
		esetName1<-paste(data1,"_",platform1,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
		if(esetpath1 != "")  # esetpath present. directly load the eset from the directory
        	{
                	eset <- try(load(esetpath1))
                	if(length(grep("^Error.*",eset))){  # if eset not found at the specified location
                        	esetpath1=""  # query from database
                	}
                	else {
                        	eset <- as(get(eset),"ExpressionSet")
				d1pData <- getpData(exp=data1,db=platform1)
                        	if(sum(is.na(match(rownames(d1pData),colnames(Biobase::exprs(eset))))) > 0)  # not matching
                        	{
                                	colnames(Biobase::exprs(eset)) <- rownames(d1pData)
                        	}
                        	Biobase::pData(eset) <- d1pData
                        	assign(esetName1,eset)
				fil <- filter_eset(get(esetName1),filterchk1,filters1,includeORexclude1)
       				ExprData1 <- Biobase::exprs(fil)
        			d1pData <- Biobase::pData(fil)
				probes <- rownames(Biobase::exprs(eset))
				geneprobes1 <- get_gene_probe_eset(platform1,unique(probes))
				mm <- match(geneprobes1[,"probe_name"],rownames(ExprData1))
				if(sum(is.na(mm)) > 0) {
                			mm <- mm[which(is.na(mm) == FALSE)]
        			}
        			ExprData1 <- ExprData1[mm,]
				matchGenes <- match(rownames(ExprData1),geneprobes1[,"probe_name"])
				rownames(ExprData1) <- geneprobes1[matchGenes,"gene"]
                	}
        	}	
		if(esetpath1 == "")  # eset path info not provided then obtain the data by querying the database
        	{
			#d1pData<-filter(filterchk1,filters1,parsedSData=d1pData,includeORexclude1)
			#d1pData<-filter(filters1,filterchk1,parsedSData=d1pData,includeORexclude1)
			d1Values<-getDataForGRS(exp=data1,filterchk1,filters1,includeORexclude1,rownames(d1pData),platform1)
			d1Values<-d1Values[which(d1Values[,"Gene"] != "" ),]
        		NoOfProbes <- dim(d1Values)[1]/dim(d1pData)[1]
        		geneNames1 <- d1Values[1:NoOfProbes,"Gene"]

       			ExprData1<-matrix(d1Values[,"Value"],ncol=dim(d1pData)[1])
       			rownames(ExprData1) <- as.vector(geneNames1)
        		ExprData1 <- as.data.frame(ExprData1)
		}
		cat("2",file=paste("/data/srv/www/htdocs/tmp/",statusFile,sep=""))
		#pvalue1 <- calculatePvalue(d1pData,d1Values,prop1,oneTailFirst)
		pvalue1 <- calculatePvalue(d1pData,ExprData1,prop1,oneTailFirst)
		cat("3",file=paste("/data/srv/www/htdocs/tmp/",statusFile,sep=""))
		if(length(grep("^Error.*",pvalue1[1]))){
			print(paste("Error: Error in computing p-values for query dataset",pvalue1[1]))
			return(paste("Error: Error in computing p-values for query dataset",pvalue1[1]))
		}
		x <- data.frame(names(pvalue1[[1]]),pvalue1[[1]]) 
		xup <- data.frame(names(pvalue1[[2]]),pvalue1[[2]]) 
		xdn <- data.frame(names(pvalue1[[3]]),pvalue1[[3]]) 
	}
	else ## user submitted file of Pvalues
	{
		x <- try(read.table(file = paste(src1,"/",data1,sep=""),header=TRUE))
		if(!is.data.frame(x)) {
			return("Error: File format error")
		}	
		xup <- data.frame()
		xdn <- data.frame()

	}
	if(src2 == "DB") 
	{
		d2Values = NULL
                ExprData2 = NULL
                esetName2 <- paste(data2,"_",platform2,"_",paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep=""),sep="")
                if(esetpath2 != "")  # esetpath present. directly load the eset from the directory
                {
                        eset <- try(load(esetpath2))
                        if(length(grep("^Error.*",eset))){  # if eset not found at the specified location
                                esetpath2=""  # query from database
                        }
                        else {
                                eset <- as(get(eset),"ExpressionSet")
				d2pData<-getpData(exp=data2,db=platform2)
                                if(sum(is.na(match(rownames(d2pData),colnames(Biobase::exprs(eset))))) > 0)  # not matching
                                {
                                        colnames(Biobase::exprs(eset)) <- rownames(d2pData)
                                }
                                Biobase::pData(eset) <- d2pData
                                assign(esetName2,eset)
                                fil <- filter_eset(get(esetName2),filterchk2,filters2,includeORexclude2)
                                ExprData2 <- Biobase::exprs(fil)
                                d2pData <- Biobase::pData(fil)
				probes <- rownames(Biobase::exprs(eset))
                                geneprobes2 <- get_gene_probe_eset(platform2,unique(probes))
                                mm <- match(geneprobes2[,"probe_name"],rownames(ExprData2))
                                if(sum(is.na(mm)) > 0) {
                                        mm <- mm[which(is.na(mm) == FALSE)]
                                }
                                ExprData2<-ExprData2[mm,]
                                matchGenes <- match(rownames(ExprData2),geneprobes2[,"probe_name"])
                                rownames(ExprData2) <- geneprobes2[matchGenes,"gene"]

                        }
                }
                if(esetpath2 == "")  # eset path info not provided then obtain the data by querying the database
                {
                        #d2pData<-filter(filterchk2,filters2,parsedSData=d2pData,includeORexclude2)
			#d2pData<-filter(filters2,filterchk2,parsedSData=d2pData,includeORexclude2)   # correction because of flipped values being passed to the function
                        d2Values <- getDataForGRS(exp=data2,filterchk2,filters2,includeORexclude2,rownames(d2pData),platform2)
                        d2Values <- d2Values[which(d2Values[,"Gene"] != "" ),]
                        NoOfProbes <- dim(d2Values)[1]/dim(d2pData)[1]
                        geneNames2 <- d2Values[1:NoOfProbes,"Gene"]

                        ExprData2<-matrix(d2Values[,"Value"],ncol=dim(d2pData)[1])
                        rownames(ExprData2) <- as.vector(geneNames2)
                        ExprData2<-as.data.frame(ExprData2)
                }

		cat("3",file=paste("/data/srv/www/htdocs/tmp/",statusFile,sep=""))
		#pvalue2 <- calculatePvalue(d2pData,d2Values,prop2,oneTailSecond)
		pvalue2 <- calculatePvalue(d2pData,ExprData2,prop2,oneTailSecond)
		cat("4",file=paste("/data/srv/www/htdocs/tmp/",statusFile,sep=""))
		if(length(grep("^Error.*",pvalue2[1]))){
			print(paste("Error: Error in computing p-values for reference dataset",pvalue1[1]))
			return("Error: Error in computing p-values for reference dataset")
		}
		y <- data.frame(names(pvalue2[[1]]),pvalue2[[1]]) 
		yup <- data.frame(names(pvalue2[[2]]),pvalue2[[2]]) 
		ydn <- data.frame(names(pvalue2[[3]]),pvalue2[[3]]) 
		
	}
        else ## user submitted file of Pvalues
        {
		y<-try(read.table(file = paste(src2,"/",data2,sep=""),header=TRUE))
		if(!is.data.frame(y)) {
                        return("Error: File format error")
                }

		yup<-data.frame()
		ydn<-data.frame()

        }

	# contruct a table with gene lists and Pvalues
	colnames(x)<-c("GeneID","pValue")
	colnames(y)<-c("GeneID","pValue")

	if(!is.null(xup) & nrow(xup)>0) colnames(xup)<-c("GeneID","pValue")
	if(!is.null(xdn) & nrow(xdn)>0) colnames(xdn)<-c("GeneID","pValue")
	if(!is.null(yup) & nrow(yup)>0) colnames(yup)<-c("GeneID","pValue")
	if(!is.null(ydn) & nrow(ydn)>0) colnames(ydn)<-c("GeneID","pValue")


	#check for homologous genes
	hg1<-Sys.time()
	if(specie1 != specie2)
	{
		if(src1 == "DB")
			homoloGenes <- find_hgenes_GRS(exp=data1,org=specie1,glist=paste(y[,"GeneID"],collapse=",",sep = ""))
		else
			homoloGenes <- find_hgenes_GRS(exp=NULL,org=specie1,glist=paste(y[,"GeneID"],collapse=",",sep = ""))

		m <- match(y[,"GeneID"],homoloGenes[,"geneID"])
		y[,"GeneID"] <- homoloGenes[m,"homoloGeneID"]
		if(nrow(yup)>0) {
			yup[,"GeneID"] <- homoloGenes[m,"homoloGeneID"]
		}
		if(nrow(ydn)>0) {
			ydn[,"GeneID"] <- homoloGenes[m,"homoloGeneID"]
		}
		
		index <- which(!is.na(y[,"GeneID"]))
		y <- y[index,]
		yup <- yup[index,]
		ydn <- ydn[index,]

	}	
	hg2<-Sys.time()
	print(paste("Start time for hg",hg1))
	print(paste("End time hg",hg2))
	print(paste("total time for hg",hg2-hg1))


	grt1<-Sys.time()
	resultdf<- try(CLEAN::GRS(x,y,tolerateWarnings = FALSE) )
	## check for errors
        #if(length(resultdf) == 1 & resultdf == "number of common gene IDs is < 2"){
        if(length(resultdf) == 1){
		return("Error: No common genes found")
	}
	if(dim(resultdf$geneTable)[1]==0){
		return("Error: No common genes found")
	}
	Xup_Yup_resultdf<- NULL
	Xup_Ydn_resultdf <- NULL
	Xdn_Yup_resultdf <- NULL
	Xdn_Ydn_resultdf <- NULL
	Xup_Y_resultdf<- NULL
	Xdn_Y_resultdf<- NULL
	X_Yup_resultdf<- NULL
	X_Ydn_resultdf<- NULL

	# Comparison for up and down regulated genes
	if(oneTailFirst & oneTailSecond)
	{
		if(nrow(xup)>0 & nrow(yup)>0)
			Xup_Yup_resultdf<- try(CLEAN::GRS(xup,yup,tolerateWarnings = FALSE) )
		if(nrow(xup)>0 & nrow(ydn)>0)
			Xup_Ydn_resultdf<- try(CLEAN::GRS(xup,ydn,tolerateWarnings = FALSE) )
		if(nrow(xdn)>0 & nrow(yup)>0)
			Xdn_Yup_resultdf<- try(CLEAN::GRS(xdn,yup,tolerateWarnings = FALSE) )
		if(nrow(xdn)>0 & nrow(ydn)>0)
			Xdn_Ydn_resultdf<- try(CLEAN::GRS(xdn,ydn,tolerateWarnings = FALSE) )
	}

	# Comparison for query up/down vs reference
	if(oneTailFirst)
	{
		if(nrow(xup)>0)
			Xup_Y_resultdf <- try(CLEAN::GRS(xup,y,tolerateWarnings = FALSE) )
		if(nrow(xdn)>0)
			Xdn_Y_resultdf <- try(CLEAN::GRS(xdn,y,tolerateWarnings = FALSE) )

	}
	
	# Comparison for query vs reference up/down
	if(oneTailSecond)
	{
		if(nrow(yup)>0)
			X_Yup_resultdf <- try(CLEAN::GRS(x,yup,tolerateWarnings = FALSE) )
		if(nrow(ydn)>0)
			X_Ydn_resultdf <- try(CLEAN::GRS(x,ydn,tolerateWarnings = FALSE) )

	}
	
	grt2<-Sys.time()
	cat("5",file=paste("/data/srv/www/htdocs/tmp/",statusFile,sep=""))

	print(paste("Start time for grs",grt1))
	print(paste("End time grs",grt2))
	print(paste("total time for grs",grt2-grt1))
	#####Generating sessionID
	sessionID<-paste(unlist(strsplit(gsub(":"," ",date())," ") ),collapse="_",sep="")
	resultdf$sessionID<-sessionID
	
	GRSName<-paste("GRSresultdf",sessionID,sep="")
	
	assign(GRSName,list(resultdf,Xup_Yup_resultdf,Xup_Ydn_resultdf,Xdn_Yup_resultdf,Xdn_Ydn_resultdf,Xup_Y_resultdf,Xdn_Y_resultdf,X_Yup_resultdf,X_Ydn_resultdf,levels1,levels2,prop1,prop2))
	save(list=GRSName,file=paste("/data/srv/www/htdocs/tmp/",GRSName,".RData",sep=""))

	x<-lapply(get(GRSName)[1:9],function(x) x$geneTable)
	allGenes<-rownames(x[[1]])
	newgeneTable<-allGenes
	for(i in 1:length(x))
	{
		if(is.null(x[[i]]))
		{
			colg1<-rep(NA,length(allGenes))
		}
		else
		{
			g1<-rownames(x[[i]])
			m<-match(allGenes,g1)
			colg1<-x[[i]][,2][m]
		}
		if(is.null(newgeneTable))
			newgeneTable<-colg1
		else
			newgeneTable<-data.frame(newgeneTable,colg1)
	}
	colnames(newgeneTable) <- c("Gene","Query_Reference","Queryup_ReferenceUp","Queryup_ReferenceDn","QueryDn_ReferenceUp","QueryDn_ReferenceDn","Queryup_Reference","QueryDn_Reference","Query_ReferenceUp","Query_ReferenceDn")
	write.table(newgeneTable, file=paste("/data/srv/www/htdocs/tmp/GRSgenetable",sessionID,".xls",sep=""),row.names=FALSE, sep="\t", quote=TRUE)

	print("processing done")
	g<-Sys.time()
	print(paste("Start time for master",f))
	print(paste("End time master",g))
	print(paste("total time for master",g-f))
	
	# print result to logs

	print(paste("X_Y_resultdf",if(!is.null(resultdf)) resultdf$p.value) )
	print(paste("Xup_Yup_resultdf",if(!is.null(Xup_Yup_resultdf)) Xup_Yup_resultdf$p.value) )
	print(paste("Xup_Ydn_resultdf",if(!is.null(Xup_Ydn_resultdf)) Xup_Ydn_resultdf$p.value) )
	print(paste("Xdn_Yup_resultdf",if(!is.null(Xdn_Yup_resultdf)) Xdn_Yup_resultdf$p.value) )
	print(paste("Xdn_Ydn_resultdf",if(!is.null(Xdn_Ydn_resultdf)) Xdn_Ydn_resultdf$p.value) )
	print(paste("X_Yup_resultdf",if(!is.null(X_Yup_resultdf)) X_Yup_resultdf$p.value) )
	print(paste("X_Ydn_resultdf",if(!is.null(X_Ydn_resultdf)) X_Ydn_resultdf$p.value) )
	print(paste("Xup_Y_resultdf",if(!is.null(Xup_Y_resultdf)) Xup_Y_resultdf$p.value) )
	print(paste("Xdn_Y_resultdf",if(!is.null(Xdn_Y_resultdf)) Xdn_Y_resultdf$p.value) )
	return(get(GRSName))
}
