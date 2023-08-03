
#' A function for ...
#'
#' This function allows you to ...
#' @param basement A two value vector containing "baseline" and "treatment".
#' @param authors M. Fazel-Najafabadi
#' @keywords test
#' @export 
#' @examples
#' res <- batchTregGRS(....)

createSignatureByDiffGEA <- function(exp, prop=NULL, filterchk="n", basement=NULL, includeORexclude=2, ifCluster="Genes", up=4000, down=1000, window.size=50,
				   path_to_write, pORqVal, userUploadedProfilePath=NULL, previousID=NULL, reAnalyse=FALSE, 
				   glist=NULL, pValue=NULL, foldchange=NULL, MLinfo=NULL, write=FALSE, homoloGenes=TRUE, multigroup=FALSE,
				   sigDB="ilincs_sigs", chost="ilincs.org", debuging=FALSE)
				   
{
#print(paste("exp",1,exp))
#print(paste("prop",1,prop))
#print(paste("filterchk",1,filterchk))
#print(paste("basement",1,basement))
#print(paste("includeORexclude",1,includeORexclude))
#print(paste("ifCluster",1,ifCluster))
#print(paste("path_to_write",1,path_to_write))
#print(paste("pORqVal",1,pORqVal))
#print(paste("userUploadedProfilePath",1,userUploadedProfilePath))
#print(paste("reAnalyse",1,reAnalyse))
#print(paste("previousID",1,previousID))
#print(paste(MLinfo))
## un-necessary part for new ilincs portal muted by #@@

	if(is.null(exp)) {
		message("missing exp!!!")
		res <- getDefaultErrorDataframe("!!! .... exp should not be missing/NULL")
		return(res)
	} else if (!is.null(userUploadedProfilePath)) {
		message("missing library!!!")
		res <- getDefaultErrorDataframe("!!! .... libName should not be missing/NULL")
		return(res)
	
	} else {
#@@		if(!is.null(userUploadedProfilePath)){
			#upload custom profile from a file
#@@			queryProfile<-try(read.table(file = paste(userUploadedProfilePath,"/",exp,sep=""),header=TRUE))
#@@			if(dim(queryProfile)[1] == 0){
#@@				res <- getDefaultErrorDataframe("empty file")
#@@				return(res)
#@@			}
			
#@@			if(dim(queryProfile)[2] == 4) colnames(queryProfile) <- c("geneID","geneName","coefficients","Pvals")
#@@			if(dim(queryProfile)[2] == 4) queryProfile <- queryProfile[, c(1,3,4)]
#@@			if(dim(queryProfile)[2] == 3) colnames(queryProfile) <- c("geneID","coefficients","Pvals")
#@@			if(dim(queryProfile)[2] == 2) colnames(queryProfile) <- c("geneID","coefficients")
#@@			if((dim(queryProfile)[2] %in% c(2,3)) & is.na(as.numeric(queryProfile[1,1]))) {
#@@			    queryProfile$geneID <- symbol2geneid(paste(queryProfile$geneID, collapse=","))$GeneID
#@@			    }
#@@			if(dim(queryProfile)[2] != 2 & dim(queryProfile)[2] != 3 & dim(queryProfile)[2] != 4){
#@@				res <- getDefaultErrorDataframe("wrong format")
#@@				return(res)
#@@			}
#@@			queryProfile <- na.omit(queryProfile)
			#####Generating sessionID
#@@		        sessionID <- paste(c(strsplit(date(),split=" ")[[1]],as.integer(runif(1)*10e6)),sep="",collapse="_")
#@@			sessionID <- gsub(":","-", sessionID)
#@@			sessionID <- gsub("__","_", sessionID)

#@@			res <- getDefaultErrorDataframe("Done")  #dummy object
#@@			res[1] <- sessionID
#@@			res <- res[-(c(16,17,18))]
#@@		} else{
	        	##@Mehdi Generating diff Gene list/profile 
	        	res <- master_diffGenes(exp=exp, prop=prop, filterchk=filterchk, basement=basement, includeORexclude=includeORexclude, 
						ifCluster=ifCluster, up=4000, down=1000, window.size=50, 
						path_to_write=path_to_write, pORqVal=pORqVal, previousID=previousID, reAnalyse=reAnalyse, 
						glist=glist, pValue=pValue, foldchange=foldchange, MLinfo=MLinfo, homoloGenes=homoloGenes, multigroup=multigroup,
						sigDB=sigDB, chost=chost, debuging=debuging)

			if(!multigroup) {
			    level1str <- res[[13]]
			    print(str(res))
			    if(as.character(res[[2]])== "NA"){ #error in anova or no genes found etc.
				res <- getDefaultErrorDataframe(res[[6]])  # remark from function
				res[[13]] <- level1str    #different if "no two levels"
				print(res)
				return(res)
			    }
			}
#@@			queryProfile <- data.frame(geneID=unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[13])))," ")), 
			#geneName=unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[14])))," ")), 
#@@			coefficients=as.numeric(unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[12])))," "))), 
#@@			Pvals=as.numeric(unlist(strsplit(unlist(gsub("DELIM"," ",unlist(res[17])))," "))))

#@@		}
		#Merge multiple rows for the same geneIds
#@@		if(dim(queryProfile)[2] == 3) unqIdPvals <- split(queryProfile[,"Pvals"], queryProfile[,"geneID"])
#@@		if(dim(queryProfile)[2] == 3) mergePvals <- unlist(lapply(unqIdPvals, function(x) exp(mean(log(x), na.rm = TRUE))))
#@@ 		unqIdCoeffs <- split(queryProfile[,"coefficients"], queryProfile[,"geneID"])
#@@		mergeCoeffs <- unlist(lapply(unqIdCoeffs, function(x) mean(as.numeric(x), na.rm = TRUE)))
#@@		if(dim(queryProfile)[2] == 3) queryProfile <- data.frame(geneID=names(unqIdCoeffs),coefficients=mergeCoeffs,Pvals=mergePvals)
#@@                if(dim(queryProfile)[2] == 2) queryProfile <- data.frame(geneID=names(unqIdCoeffs),coefficients=mergeCoeffs)
## Mehdi: have to turn it on if need geneNames, this way is slow, we neet to transfer them not query
#                queryProfile$geneName <- geneid2symbol(paste(queryProfile$geneID, collapse=","))$Symbol
#                queryProfile <- queryProfile[,c("geneID", "geneName", "coefficients", "Pvals")]
#@@		write.table(queryProfile, file=paste(path_to_write,"queryProfile_",as.character(res[[1]]),".xls",sep=""), 
#@@			    sep="\t", row.names=FALSE) #, row.names=F, col.names=F,append=T)
	}
	#### Start Parsing 'res' to test Volcano plot
# # # # #     if (!multigroup) {
# # # # # 		parse_res(res, path_to_write, write=write)
# # # # #     }

    if (reAnalyse) {
      print(paste("This action is removed...!"))
# # # # #       file.copy(paste(path_to_write,"temp_volcano_", previousID, ".xls", sep = ""), 
# # # # # 		paste(path_to_write,"temp_volcano_", as.character(res[[1]]), ".xls", sep = ""), overwrite = TRUE)
# # # # #       file.copy(paste(path_to_write,"temp_volcano_", previousID, ".RData", sep = ""), 
# # # # # 		paste(path_to_write,"temp_volcano_", as.character(res[[1]]), ".RData", sep = ""), overwrite = TRUE)
    }
    
	'if (res$coefficientsGeneID != "NA" & res$coefficients != "NA" & res$Pvals != "NA") {# & res$gpnames != "NA") {
        
        GeneID <- strsplit(res$coefficientsGeneID, split = "DELIM")
        GeneID <- as.vector(GeneID[[1]])

        GeneNames <- strsplit(res$coefficientsGeneName, split = "DELIM")
        GeneNames <- as.vector(GeneNames[[1]])
        
        coefficients <- strsplit(res$coefficients, split = "DELIM")
        coefficients <- as.vector(coefficients[[1]])
        coefficients <- as.numeric(coefficients)
        
        Pvals <- strsplit(res$Pvals, split = "DELIM")
        Pvals <- as.vector(Pvals[[1]])
        Pvals <- as.numeric(Pvals)
        
        SignatureTable <- data.frame(GeneID, GeneNames, coefficients, Pvals, stringsAsFactors = F)
        volcano_name <- paste(path_to_write,"temp_volcano_", as.character(res[[1]]), ".RData", sep = "")
        volcano_name2 <- paste(path_to_write,"temp_volcano_", as.character(res[[1]]), ".txt", sep = "")
        save(SignatureTable, file = volcano_name)
        write.table(SignatureTable, file = volcano_name2, row.names=FALSE, sep="\t")
        }
      }
      '
        #### End parsing 'res' to test Volcano plot

	return(res) 
}
