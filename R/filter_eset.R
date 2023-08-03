
########## This is the new function to filter any "eset" based on choosen property #############
## eset shoul be in class of ExpressionSet (Biobase) 					      ##
## filterchk is a string with no spaces, property are separated by ":" and values are 	      ##
## separated by ",". example: filterchk <- "subtype:Basal,subtype:Luminal,ER:+,ER:NA"	      ##
## if the path_to_write is provided then an .RData file containing the eset with the specified##
## name will be saved. If name is not provided, the default name is "filtered.eset".	      ##
## Written by Mehdi.Fazel								      ##
################################################################################################

#' A function for filtering Biobase ExpressionSets.
#'
#' This function allows you to filter an ExpressionSet based on specified properties.
#' @param eset should be an ExpressionSet with the class of Biobase Eset containing pData and exprs data
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
#' @param new.esetName This parameter should be the eset name, if user wants to have a specific name.
#' @param path_to_write This optional parameter specifies the path which user wants to save the 
#'	filtered expressionSet to be saved in as an .RData file.
#' @param authors M. Fazel-Najafabadi
#' @keywords filter
#' @export 
#' @examples
#' ## not run
#' data(expl_eset)
#' res <- filter_eset(eset=expl_eset, filterchk="ER:+,ER:-", includeORexclude=1)
#' ## not run


filter_eset <- function(eset, filterchk=NULL, basement=NULL, includeORexclude,
			new.esetName="filtered.eset", path_to_write=NULL) {
#library(Biobase)
if (class(eset)[1] == "exprSet") eset <- as(eset, "ExpressionSet")
if (is.null(filterchk)) return(assign(new.esetName, eset))
if (filterchk == "n") return(assign(new.esetName, eset))
filterchk <- gsub("<NA>", "NA", filterchk)

pdat <- as.data.frame(Biobase::pData((eset)), useNA="ifany")#, stringsAsFactors=FALSE)
# # # for(j in colnames(pdat)) pdat[,j] <- as.character(pdat[,j])
# # # pdat <- data.frame(apply(pdat, 2, as.character), stringsAsFactors=F)
##@i <- sapply(pdat, is.factor)
# pdat[i] <- lapply(pdat[i], as.character)
##@if(sum(i) > 0) pdat[,which(i==T)] <- lapply(pdat[,which(i==T)], as.character)
pdat[is.na(pdat)] <- "NA"
pdat$MeasurementName <- rownames(pdat)
edat <- Biobase::exprs((eset)) ; fdat <- Biobase::fData((eset)) # ; fdat <- as.data.frame(apply(fdat, 2, as.character), stringsAsFactors=FALSE)

    prop <- as.data.frame(t(sapply(unlist(strsplit(filterchk, ",,,")), function(i) unlist(strsplit(i, ":")))), 
			  row.names=NULL, stringsAsFactors=FALSE)
    if (sum(unique(prop[,1]) %in% colnames(pdat)) < length(unique(prop[,1]))) return(message("Filtering property out of range...!"))
	if(as.character(includeORexclude) == "1"){ ## this stays!!!
    #	    for (pro in unique(prop[,1])) pdat <- pdat[which(pdat[, pro] %in% as.character(prop[which(prop[,1]==pro),2])),]
		for (pro in unique(prop[,1])) pdat <- pdat[(pdat[, pro] %in% prop[which(prop[,1]==pro),2]), ,drop=F]
		edat <- as.matrix(edat[, rownames(pdat), drop=FALSE])
    ##		pdat <- pdat[lastpro, ]; edat <- as.matrix(edat[, lastpro])
	    } else { ## this stays!!!
    #	    for (pro in unique(prop[,1])) pdat <- pdat[-which(pdat[, pro] %in% (prop[which(prop[,1]==pro),2])),]
		for (pro in unique(prop[,1])) pdat <- pdat[!(pdat[, pro] %in% prop[which(prop[,1]==pro),2]), ,drop=F]
		    edat <- as.matrix(edat[, rownames(pdat), drop=FALSE])
}
assign(new.esetName, new("ExpressionSet", exprs=as.matrix(edat), 
	  phenoData=new("AnnotatedDataFrame", data=pdat), 
	  featureData=new("AnnotatedDataFrame", data=as.data.frame(fdat))))
if (!is.null(path_to_write)) {
	  save(list=new.esetName, file=paste(path_to_write, new.esetName, ".RData",sep=""))}
return(get(new.esetName))
}
