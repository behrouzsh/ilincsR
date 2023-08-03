
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

create_jnlps <- function(sessionID, isCentered, path_to_write, sigDB="ilincs_sigs", chost="ilincs.org", debuging=FALSE){ #, ilincs=TRUE){

jnlpF <- c(
	    "<?xml version=\"1.0\" encoding=\"utf-8\"?>", 							# 1
	    "<!-- JNLP File for Notepad -->", 									# 2
	    "<jnlp spec=\"1.0+\"", 										# 3
		"\tcodebase=\"http://eh3.uc.edu/clean/fTreeView\"", 						# 4
		"\tinsert", 											# 5
		"\t<information>", 										# 6
		    "\t\t<title>Enhanced Java TreeView</title>", 						# 7
		    "\t\t<vendor>", 										# 8
		    "\t\t\tLaboratory for Statistical Genomics, ", 						# 9
		    "\t\t\tUniv. Cincinnati.", 									# 10
		    "\t\t</vendor>", 										# 11
		    "\t\t<description kind=\"one-line\">", 							# 12
			"\t\t\tJava TreeView with enhanced features.", 						# 13
		    "\t\t</description>", 									# 14
		    "\t\t<description kind=\"tooltip\">", 							# 15
		    "\t\t\tEnhanced Java TreeView", 								# 16
		    "\t\t</description>", 									# 17
		    "\t\t<homepage href=\"http://eh3.uc.edu\"/>", 						# 18
		    "\t\t<offline-allowed/>", 									# 19
# 		    "\t\t<!--", 										# 20
# 		    "\t\t<shortcut online=\"false\">", 								# 21
# 			"\t\t\t<desktop/>", 									# 22
# 		    "\t\t\t<menu submenu=\"TreeView\"/>", 							# 23
# 		    "\t\t</shortcut>", 										# 24
# 		    "\t\t-->", 											# 25
		"\t</information>", 										# 26
		"\t<update check=\"always\" policy=\"prompt-update\"/>", 					# 27
		"\t<security>", 										# 28
		    "\t\t<all-permissions/>", 									# 29
		"\t</security>",  										# 30
		"\t<resources>\t", 										# 31
		    "\t\t<j2se version=\"1.5+\"/>", 								# 32
		    "\t\t<jar href=\"jars/LFTreeView.jar\" main=\"true\"/>", 					# 33
# 		    "\t\t<!--", 										# 34
# 		    "\t\t<jar href=\"jars/plugins/Dendrogram.jar\"/>", 						# 35
# 		    "\t\t<jar href=\"jars/nanoxml-2.2.2.jar\"/>", 						# 36
# 		    "\t\t-->", 											# 37
		"\t</resources>", 										# 38
		"\t<application-desc main-class=\"edu.stanford.genetics.treeview.app.LinkedViewApp\">",  	# 39
		    "\t\t<argument>-r</argument>", 								# 40
		    "\t\tinsert", 										# 41
		    "\t\t<argument>-new</argument>", 								# 42
		"\t</application-desc>", 									# 43
	    "</jnlp>")												# 44

# jnlpF <- c("<?xml version=\"1.0\" encoding=\"utf-8\"?>", "<!-- JNLP File for Notepad -->", "<jnlp spec=\"1.0+\"",
#  "\tcodebase=\"http://eh3.uc.edu/clean/fTreeView/\"", "\tinsert", "\t<information>", "\t\t<title>Enhanced Java TreeView</title>",
#  "\t\t<vendor>", "\t\t\tLaboratory for Statistical Genomics, ", "\t\t\tUniv. Cincinnati.", "\t\t</vendor>",
#  "\t\t<description kind=\"one-line\">", "\t\t\tJava TreeView with enhanced features.", "\t\t</description>",
#  "\t\t<description kind=\"tooltip\">", "\t\t\tEnhanced Java TreeView", "\t\t</description>",
#  "\t\t<homepage href=\"http://eh3.uc.edu\"/>", "\t\t<offline-allowed/>", "\t\t<!--", "\t\t<shortcut online=\"false\">",
#  "\t\t\t<desktop/>", "\t\t\t<menu submenu=\"TreeView\"/>", "\t\t</shortcut>", "\t\t-->", "\t</information>",
#  "\t<update check=\"always\" policy=\"prompt-update\"/>", "\t<security>", "\t\t<all-permissions/>", "\t</security>",
#  "\t<resources>\t", "\t\t<j2se version=\"1.5+\"/>", "\t\t<jar href=\"LFTreeView.jar\" main=\"true\"/>",
#  "\t\t<!--", "\t\t<jar href=\"jars/plugins/Dendrogram.jar\"/>", "\t\t<jar href=\"jars/nanoxml-2.2.2.jar\"/>",
#  "\t\t-->", "\t</resources>", "\t<application-desc main-class=\"edu.stanford.genetics.treeview.app.LinkedViewApp\">",
#  "\t\t<argument>-r</argument>", "\t\tinsert", "\t\t<argument>-new</argument>", "\t</application-desc>", "</jnlp>")

#jnlpF <- scan("http://www.eh3.uc.edu/tmp/TreeView.jnlp.temp",what="character",sep="\n",quote="",allowEscapes=F)
	
# 	machine <- system("hostname", intern=TRUE)
# 	if(machine == "GIMM14") loc <- "http://eh3.uc.edu" else loc <- "http://dev.ilincs.org"
	servSet <- getServerSettings(debuging=debuging, sigDB=sigDB, chost=chost)
	
	if(servSet$compute == "ilincs.org") {
	    loc <- "http://eh3.uc.edu" 
	} else if(servSet$compute == "dev.ilincs.org") {
	    loc <- "http://dev.ilincs.org"
	} else return("Unknown server!!!")
	if(isCentered) {
	    jnlpFNC <- paste0(path_to_write,"TreeView-Centered", sessionID, ".jnlp")
	    cat(jnlpF[1:4], file=jnlpFNC, sep="\n")
# 	    jnlpLine <- paste0("\thref=\"http://eh3.uc.edu/tmp/TreeView-Centered", sessionID, ".jnlp\">")
	    jnlpLine <- paste0("\thref=\"", loc, "/tmp/TreeView-Centered", sessionID, ".jnlp\">")
	    cat(jnlpLine, file=jnlpFNC,sep="\n", append=T)
	    cat(jnlpF[6:30], file=jnlpFNC,sep="\n", append=T) #6:40
	    #jnlpLine<-paste("<argument>","http://eh3.uc.edu/tmp/TreeView-Centered",sessionID,".cdt","<\/argument>",sep="")
# 	    jnlpLine <- paste(paste0("\t", "<argument>","http://eh3.uc.edu/tmp/TreeView-Centered", sessionID, ".cdt","</argument>"), 
	    jnlpLine <- paste(paste0("\t", "<argument>", loc, "/tmp/TreeView-Centered", sessionID, ".cdt","</argument>"), 
## 	    jnlpLine <- paste(paste("\t", "<argument>", path_to_write, "TreeView-Centered",sessionID,".cdt","</argument>",sep=""), 
	    paste0("\t", "<argument>-overwrite_global_config</argument>"), sep="\n")
	    cat(jnlpLine, file=jnlpFNC, sep="\n", append=T)
# 	    jnlpLine <- paste("<argument>-overwrite_global_config</argument> ")
# 	    cat(jnlpLine,file=jnlpFNC,sep="\n",append=T)
	    cat(jnlpF[32:34], file=jnlpFNC, sep="\n", append=T) #42:44
# 	    if(ilincs) cat("\t","\t","<argument>-overwrite_global_config</argument>",file=jnlpFN,sep="\n",append=T)
# 	    if(ilincs) cat("\t","\t","<argument>-ilincs</argument>","\n",file=jnlpFN,append=T)
	
	} else {
	    jnlpFN <- paste0(path_to_write, "TreeView", sessionID,".jnlp")
	    cat(jnlpF[1:4], file=jnlpFN, sep="\n")
# 	    jnlpLine <- paste0("\thref=\"http://eh3.uc.edu/tmp/TreeView", sessionID, ".jnlp\">")
	    jnlpLine <- paste0("\thref=\"", loc, "/tmp/TreeView", sessionID, ".jnlp\">")
	    cat(jnlpLine, file=jnlpFN, sep="\n", append=T)
	    cat(jnlpF[6:30], file=jnlpFN, sep="\n", append=T) #6:40
# 	    jnlpLine<-paste("<argument>","http://eh3.uc.edu/tmp/TreeView",sessionID,".cdt","<\/argument>",sep="")
# 	    jnlpLine <- paste(paste0("\t", "<argument>","http://eh3.uc.edu/tmp/TreeView", sessionID, ".cdt","</argument>"), 
	    jnlpLine <- paste(paste0("\t", "<argument>", loc, "/tmp/TreeView", sessionID, ".cdt","</argument>"), 
##	    jnlpLine <- paste(paste("\t", "<argument>", path_to_write, "TreeView",sessionID,".cdt","</argument>",sep=""), 
	    paste0("\t", "<argument>-overwrite_global_config</argument>"), sep="\n")
	    cat(jnlpLine, file=jnlpFN, sep="\n", append=T)
# 	    jnlpLine <- paste("<argument>-overwrite_global_config</argument> ")
# 	    cat(jnlpLine,file=jnlpFNC,sep="\n",append=T)
	    cat(jnlpF[32:34], file=jnlpFN, sep="\n", append=T) #42:44
# 	    if(ilincs) cat("\t","\t","<argument>-overwrite_global_config</argument>",file=jnlpFN,sep="\n",append=T)
# 	    if(ilincs) cat("\t","\t","<argument>-ilincs</argument>","\n",file=jnlpFN,append=T)
	}
# print(paste(jnlpFN))
}

