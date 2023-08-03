
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

getProbes <- function(geneid=NULL,experiment=NULL,db=NULL, debuging=FALSE){
        probes <- NULL
	servSet <- getServerSettings(debuging)
	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = db, port = servSet$port, host = servSet$host, password = "public")
#        if (!test) {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="10.165.4.231", port = 4040, password='public')}
#	else {
#	  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname=db, host="gimm2.ketl.uc.edu", password='public')}
	#mycon <- DBI::dbConnect(RMySQL::MySQL(), user="public", dbname=platform, host="db2.ketl.uc.edu", port = 3306, password="public")

        F <- "select AT.Name from ArrayType AT, Array A,Measurement M, MeasurementList ML,Experiment E where E.Measurement_List_ID=ML.ID and ML.Measurement_ID=M.ID and M.Name= A.Name and A.Array_Type_ID=AT.ID and E.Name='"
        F1 <- "'"
        Fsql <- paste(F,experiment,F1,sep="")
        rs <- DBI::dbSendQuery(mycon, Fsql)
        AT <- DBI::fetch(rs, n = -1)

        geneid <- paste(geneid,collapse=",",sep="")
	geneid <- gsub(",,",",",geneid)
        sql <- paste("SELECT F.Name, F.ID, G.Name FROM Gene G, GeneList GL, Feature F, Reporter R WHERE G.Name IN (",geneid,")  AND G.Name NOT LIKE \"%/%\" AND GL.Gene_ID = G.ID AND R.Gene_List_ID = GL.ID AND F.Reporter_ID = R.ID",sep="")
        rs <- DBI::dbSendQuery(mycon, sql)
        probes <- DBI::fetch(rs,n=-1)
        colnames(probes) <- c("Feature_Name","Feature_ID","Gene")
        DBI::dbDisconnect(mycon)
        return(probes)
}

# getProbes <- function(geneid=NULL,experiment=NULL,db=NULL, debuging=FALSE){
#         probes <- NULL
# 	servSet <- getServerSettings(test) 
# # 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "GPDatasets", port = servSet$port, host = servSet$host, password = "public")
# 	mycon <- DBI::dbConnect(RMySQL::MySQL(), user = "public", dbname = "QueryTables", port = servSet$port, host = servSet$host, password = "public")
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
# # 	pn <- DBI::dbGetQuery(mycon, paste0("SELECT platformDb FROM ExperimentMetadata WHERE experimentName = '", experiment, "'"))[1,1]
# 	pn <- DBI::dbGetQuery(mycon, paste0("SELECT Platform FROM ExperimentMetadata WHERE Experiment = '", experiment, "'"))[1,1]
# # 	pi <- DBI::dbGetQuery(mycon, paste0("SELECT ID FROM PlatformInfo WHERE Platform = '", pn, "'"))[1,1]
# 	pi <- DBI::dbGetQuery(mycon, paste0("SELECT ID FROM Platforms WHERE Platform = '", pn, "'"))[1,1]
# # 	probes <- DBI::dbGetQuery(mycon, paste0("SELECT Probe, Feature_ID, Gene FROM ProbeGene WHERE Platform_ID = '", pi, "' AND Gene IN (",geneid,")"))
# 	probes <- DBI::dbGetQuery(mycon, paste0("SELECT Probe, Gene FROM ProbeGene WHERE Platform_ID = '", pi, "' AND Gene IN (",geneid,")"))
#         colnames(probes) <- c("Feature_Name","Gene")
#         DBI::dbDisconnect(mycon)
#         return(probes)
# }

