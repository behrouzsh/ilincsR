
#' A function for ...
#'
#' This function allows you to ...
#' @param exp This is basically the experiment name.
#' @param test This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if test="TRUE" to make sure the uploaded dataset/experiment is working fine.
#' @param authors M. Fazel-Najafabadi
#' @keywords filter
#' @export
#' @examples
#' res <- get_metadata(....)

get_metadata <- function(exp, debuging=FALSE) {
  #exp="LDG-1193"
#  library(RMySQL)
servset <- getServerSettings(debuging=debuging)
  mycon <- DBI::dbConnect(RMySQL::MySQL(), user='public', dbname="QueryTables", port = servset$port, host=servset$host, password='public')
  #mycon <- dbConnect(MySQL(), user='public', dbname="QueryTables", port = 3306,host="gimm2.ketl.uc.edu", password='public')
#  sql <- paste("select Organism,Platform,DataFormat from ExperimentMetadata where Experiment = '",exp,"'",sep="")
  sql <- paste("select Platform from ExperimentMetadata where Experiment = '",exp,"'",sep="")
#  rs <- DBI::dbSendQuery(mycon, sql)
#  tmp <- DBI::fetch(rs,n=-1)
  #dbDisconnect(mycon)

#  dataFormat <- tmp[1,"DataFormat"]
#  org <- tmp[1,"Organism"]
#  platform <- tmp[1,"Platform"]
  platform <- DBI::dbGetQuery(mycon, sql)[1,"Platform"]
  #mycon <- dbConnect(MySQL(), user='public', dbname="QueryTables", port = 4040,host="10.165.4.231", password='public')
  #mycon <- dbConnect(MySQL(), user='public', dbname="QueryTables", port = 3306,host="gimm2.ketl.uc.edu", password='public')
  sql <- paste("select A.Value as Value from ",platform,".Experiment E,",platform,".AttributeText A where E.Attributes_ID=A.ID 
		    and E.Name='",exp,"'and The_Index=0",sep="")

#  rs <- DBI::dbSendQuery(mycon, sql)
#  tmp <- DBI::fetch(rs,n=-1)
#  DBI::dbDisconnect(mycon)
#  options <- tmp[1,"Value"]
  options <- DBI::dbGetQuery(mycon, sql)[1,"Value"]
  return(options)
}

