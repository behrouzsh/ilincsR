
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

getVariableCountsFrompData <- function(exp,db,format, debuging=FALSE)
{
##        library(DBI)
##        library(RMySQL)
        return_str <- NULL
        if(format == "Gff"){
		servSet <- getServerSettings(debuging)
		mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
#                if (!test) {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=db, port = 4040,host="10.165.4.231", password="public")}
#		else {
#		  mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=db, host="gimm2.ketl.uc.edu", password="public")}
                sql <- paste("select Distinct SP.Property,SP.Value from Dataset D,Dataset_Sample DS, Sample_Property SP where D.Dataset_Name = '",exp,"' AND DS.Dataset_ID = D.Dataset_ID AND DS.Sample_ID = SP.Sample_ID",sep = "");
                rs <- DBI::dbSendQuery(mycon, sql)
                gffprop <- DBI::fetch(rs, n=-1)
                DBI::dbDisconnect(mycon)
                parsedpData <- t(as.matrix(gffprop[,2]))
                colnames(parsedpData) <- gffprop[,1]
        } else {
                meas_limits <- initiate(db,exp,web=FALSE,dataFormat="MaxD", debuging=debuging)
                sampledata <- get_sample_properties(db,exp,meas_limits, debuging)
                parsedpData <- NULL
                parsedSD <- strsplit(sampledata[,2],";")
                #parsedSD <- data.frame(parsedSD)
                for(c in 1:length(parsedSD[[1]])){

                        tempstr <- unlist(lapply(parsedSD,function(x) x[c]))
                #       substr(str,0,regexpr('=',str)-1)
                        parsedpData <- cbind(parsedpData,substr(tempstr,regexpr('=',tempstr)+1,nchar(tempstr)))
                }
                colnames(parsedpData) <- substr(parsedSD[[1]],0,regexpr('=',parsedSD[[1]])-1)
                #parsedpData<-parse_sample_properties(propertyTags=sampledata)
        }
                variable_str <- NULL
                for(i in 1: length(colnames(parsedpData))) {
                        variable_str <- table(as.factor(parsedpData[,i]))
                        return_str <- c(return_str, paste(paste(names(variable_str),variable_str,sep= "="),collapse=","))
                }
                return_str<-paste(colnames(parsedpData),return_str,sep=":",collapse=";")
        return(return_str)
}
