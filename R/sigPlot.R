
#' A function for ...
#'
#' This function allows you to ...
#' @param sigFile This should be a signature file created or uploaded by user.
#' @param sigID This is a precomputed signature from any sig-library.
#' @param path_to_write This parameter specifies the path which user wants to save the results.
#' @param display Display plot or just get the object.
#' @param debuging This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if debuging=TRUE dev servers will be used.
#' @param authors M. Fazel-Najafabadi
#' @keywords filter
#' @export
#' @importFrom plotly %>%
#' @examples
#' ## Do not run
#' sigFile <- "sig_Thu_Apr_12_12_45_57_2018_298226.xls"
#' sigID <- "LINCSCP_109"
#' path_to_write <- "/opt/raid10/genomics/mehdi/ilincs/output/"
#' plt <- sigPlot(sigFile, sigID, path_to_write, debuging=TRUE, output=TRUE)
#' ## End do not run


sigPlot <- function(sigFile=NULL, sigID=NULL, path_to_write="/mnt/raid/tmp/", output="json", dotSize=3, fontSize=15, debuging=FALSE)
{

    res <- list(fileName=NA, Remark="Done", obj=NA)
    if(is.null(sigFile)) {
	res$Remark <- "Missing signature file!"
	return(res)
    }
    if(is.null(sigID)) {
	res$Remark <- "Missing signature ID!"
	return(res)
    }
a <- Sys.time()
print(paste("sigPlot called with following parameters at:",a))
print(paste("sigFile:", sigFile))
print(paste("sigID:", sigID))
print(paste("path_to_write:", path_to_write))
print(paste("debuging:", debuging))

x_sig <- try(read.table(file=paste0(path_to_write, sigFile), header=TRUE, stringsAsFactors=FALSE,  sep="\t")) # , quote=""

if(class(x_sig) != "data.frame") {
    res$Remark <- "No appropriate signature file or path!"
    return(res)
}

x_sig <- x_sig[, c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")]

if(dim(x_sig)[2] == 0) {
    res$Remark <- "No genes in the signature!"
    return(res)
}
genes <- paste(x_sig$ID_geneid, collapse=",")
y_sig <- downloadSignature(sigID, path_to_write=path_to_write, display=TRUE, glist=genes, geneNames=TRUE, write=FALSE, debuging=debuging)
y_sig <- y_sig[, c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")]

# dat <- merge(x=x_sig, y=y_sig, by.x=c("ID_geneid", "Name_GeneSymbol"), by.y=c("ID_geneid", "Name_GeneSymbol"))
ids <- (intersect(x_sig$ID_geneid, y_sig$ID_geneid)) # as.character

x_sig <- x_sig[x_sig$ID_geneid %in% ids, ]
dat <- data.frame(ID_geneid=x_sig$ID_geneid, Name_GeneSymbol=x_sig$Name_GeneSymbol, x_val=x_sig$Value_LogDiffExp, y_val=y_sig$Value_LogDiffExp[match(x_sig$ID_geneid, y_sig$ID_geneid)], stringsAsFactors=FALSE)

if(dim(dat)[2] == 0) return("No common genes between signatures!")
updn <- FALSE
if(abs(dat$x_val[1])==1){
  if(all(abs(dat$x_val)==1)) {
    updn <- TRUE
    dat$x_val[dat$x_val==1] <- "up"
    dat$x_val[dat$x_val==-1] <- "down"
  }
}

colnames(dat)[3:4] <- c("Created/Uploaded signature", sigID)

if("console" %in% output & "json" %in% output) {
    res$Remark <- "Can not return json and console! ;-)"
    return(res)
}

# pdf(paste0(path_to_write,"mytest.pdf"))
# plotly::plot_ly(dat, x = dat$Value_LogDiffExp.x, y = dat$Value_LogDiffExp.y, 
# 			text = dat$Name_GeneSymbol, mode = "markers") #, type="scatter") # , text = paste("Clarity: ", clarity), color = carat, size = carat)
# 
# dev.off()

# plt <- plotly::plot_ly(dat, x = ~ `Created/Uploaded signature`, y = as.formula(paste("~", names(dat)[4])), 
# 			    text = dat$Name_GeneSymbol, mode = "markers", type = "scatter")
if("png" %in% output | "pdf" %in% output | "ggplot" %in% output | "console" %in% output) {
    if(!updn) {
	plt <- ggplot2::ggplot(dat, ggplot2::aes_(x = as.name(names(dat)[3]), y = as.name(names(dat)[4])), fill=as.name(colnames(dat)[2])) + 
		ggplot2::geom_point(color = "blue", size = dotSize) + 
		ggplot2::theme_bw() +
		ggplot2::theme(text = ggplot2::element_text(size = fontSize), axis.text = ggplot2::element_text(size = fontSize))
    } else {
	plt <- ggplot2::ggplot(dat, ggplot2::aes_(x = as.name(names(dat)[3]), y = as.name(names(dat)[4])), fill=as.name(colnames(dat)[2])) + 
		ggplot2::geom_boxplot(outlier.shape = NA) + 
		ggplot2::geom_jitter(shape=16, size = dotSize, position = ggplot2::position_jitter(0.2)) + 
# 		ggplot2::geom_jitter(width=0.2, size = dotSize) + #shape=16, size = dotSize, position = ggplot2::position_jitter(0.2)) + 
		ggplot2::theme_bw() +
		ggplot2::theme(text = ggplot2::element_text(size = fontSize), axis.text = ggplot2::element_text(size = fontSize))
    }
# # # 	    ggplot2::theme_bw() +
# # # 	    ggplot2::theme(text = ggplot2::element_text(size = fontSize), axis.text = ggplot2::element_text(size = fontSize))
# # # 	    if(updn) plt <- plt + ggplot2::geom_boxplot() + ggplot2::geom_jitter(shape=16, size = dotSize, position = ggplot2::position_jitter(0.2))
#     res$obj <- plt
}
    
    sessionID <- generateSID()
    res$fileName <- paste0("sigPlot_", sessionID)

    if("png" %in% output) {
	file <- paste0(path_to_write, "sigPlot_", sessionID, ".png")
	ggplot2::ggsave(file, plot=plt, dpi=600)
    #     try(dev.off())
# # # 	return(paste0("Done file: ", path_to_write, "gseaPlot_", sessionID, ".png"))
    }
    if("pdf" %in% output) {
	file <- paste0(path_to_write, "sigPlot_", sessionID, ".pdf")
	ggplot2::ggsave(file, plot=plt, dpi=600)
    #     try(dev.off())
# # # 	return(paste0("Done file: ", path_to_write, "gseaPlot_", sessionID, ".pdf"))
    }
    if("ggplot" %in% output) {
# 	res$fileName <- plt
	save(plt, file = paste0(path_to_write, "sigPlot_", sessionID, ".RData"))
# # # 	return(res)
    }
    
    if("json" %in% output | "widget" %in% output) {
    if(!updn){
      require("plotly", quietly=TRUE)
  #     require("plotly")
  # # # 	pl <- plotly::ggplotly(plt)
	  pl <- (plotly::plot_ly(dat, x = ~ `Created/Uploaded signature`, y = as.formula(paste("~", names(dat)[4])), 
			      text = dat$Name_GeneSymbol, 
			      mode = "markers",
  # 			    jitter = ifelse(updn, 0.5, 0),
  # 			    size = dotSize,
			      marker = list(size = 10),
			      type = ifelse(updn, "box", "scatter")) %>%
	  plotly::layout(xaxis = list(title="Created/Uploaded signature", tickfont = list(size = fontSize), titlefont = list(size = fontSize)), 
			yaxis = list(title=names(dat)[4], tickfont = list(size = (fontSize+3)), titlefont = list(size = (fontSize+3)))
			))
      } else {
	  colnames(dat)[3] <- "X"
	  pl <- ggplot2::ggplot(dat, ggplot2::aes_(x=as.name(colnames(dat)[3]), y=as.name(colnames(dat)[4]), label=as.name(colnames(dat)[2]), fill=as.name(colnames(dat)[3]))) + 
		  ggplot2::geom_boxplot(outlier.shape = NA) + 
		  ggplot2::ylab("logDiffExp") +
      # 	    ggplot2::geom_point() + 
# 		  ggplot2::geom_jitter(shape=16, position = ggplot2::position_jitter(0.2)) +
		  ggplot2::geom_jitter(width=0.2) +
		  ggplot2::theme_bw() +
		  ggplot2::theme(text = ggplot2::element_text(size = fontSize), axis.text = ggplot2::element_text(size = fontSize))
	  pl <- plotly::ggplotly(pl)

      }
    }
    
    if("widget" %in% output) {
	htmlwidgets::saveWidget(pl, paste0(path_to_write, "sigPlot_", sessionID, ".html"), selfcontained = TRUE)
# 	res$obj <- pl
# # # 	return(res)
    }
    
    if("console" %in% output) {
	res$obj <- plt
    }
    if("json" %in% output) {
# # # 	pl <- plotly::ggplotly(plt)
	pl <- plotly::plotly_json(pl, FALSE)
	res$obj <- pl # plotly::plotly_json(plt, FALSE)
# 	return(res)
    }
    if(sum(output %in% c("png", "pdf", "widget", "json", "ggplot")) == 0) {
	res$Remark <- "Error!!!"
# 	return(res)
    }
    return(res)

}
