
#' A function for ...
#'
#' This function allows you to ...
#' @param sigfile ...
#' @param sigID ...
#' @param sigType ...
#' @param xlimit ...
#' @param ylimit ...
#' @param  ...
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



# # # # # volcanoPlot <- function(sigfile, sigType="notENCODE", xlimit=c(-6,6), ylimit=c(0,10), DiffExpLowerCutOff=-2, DiffExpUpperCutOff=2, 
# # # # # 			pValCutoff=0.05, labels=FALSE, cex=0.6, path_to_write="/mnt/raid/tmp/")
# # # # # {
# # # # # 
# # # # # textxy <- function (X, Y, labs, m = c(0, 0), cex = 0.5, offset = 0.8, ...)
# # # # # {
# # # # #     posXposY <- ((X >= m[1]) & ((Y >= m[2])))
# # # # #     posXnegY <- ((X >= m[1]) & ((Y < m[2])))
# # # # #     negXposY <- ((X < m[1]) & ((Y >= m[2])))
# # # # #     negXnegY <- ((X < m[1]) & ((Y < m[2])))
# # # # #     if (sum(posXposY) > 0)
# # # # #         text(X[posXposY], Y[posXposY], labs[posXposY], adj = c(0.5 - offset, 0.5 - offset), cex = cex, ...)
# # # # #     if (sum(posXnegY) > 0)
# # # # #         text(X[posXnegY], Y[posXnegY], labs[posXnegY], adj = c(0.5 - offset, 0.5 + offset), cex = cex, ...)
# # # # #     if (sum(negXposY) > 0)
# # # # #         text(X[negXposY], Y[negXposY], labs[negXposY], adj = c(0.5 + offset, 0.5 - offset), cex = cex, ...)
# # # # #     if (sum(negXnegY) > 0)
# # # # #         text(X[negXnegY], Y[negXnegY], labs[negXnegY], adj = c(0.5 + offset, 0.5 + offset), cex = cex, ...)
# # # # # }
# # # # # 
# # # # #   tbl <- read.table(paste0(path_to_write, sigfile), header=TRUE, stringsAsFactors=FALSE)
# # # # #   
# # # # #   autoScale <- TRUE
# # # # #   if (sigType=="ENCODE") {
# # # # # 
# # # # #     if (autoScale){
# # # # #       # Make a basic volcano plot
# # # # #       with(tbl, plot(Score, Probability, pch=20, main="Scatter plot")) #ENC
# # # # #     } else {
# # # # #       # Make a basic volcano plot  
# # # # #       with(tbl, plot(Score, Probability, pch=20, main="Volcano plot", xlim=xlimit,ylim=ylimit)) #ENC      
# # # # #     }
# # # # # 
# # # # # 
# # # # #     # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
# # # # #     with(subset(tbl, p.adjust(Probability, "none", length(Probability)) < pValCutoff), points(Score, Probability, pch=20, col="red"))
# # # # #     with(subset(tbl, Score >= DiffExpUpperCutOff), points(Score, Probability, pch=20, col="orange"))
# # # # #     with(subset(tbl, p.adjust(Probability,"none", length(Probability)) <= pValCutoff & (Score >= DiffExpUpperCutOff)), points(Score, Probability, pch=20, col="green"))
# # # # #     
# # # # #     # Label points with the textxy function from the calibrate plot
# # # # #     if (labels){
# # # # # #       library(calibrate)
# # # # #       with(subset(tbl, p.adjust(Probability,"none",length(Probability)) <= pValCutoff & (Score >= DiffExpUpperCutOff)), textxy(Score, Probability, labs=Name_GeneSymbol, cex=cex))
# # # # #     }
# # # # # 
# # # # #   } else {
# # # # # 
# # # # # #Example
# # # # #     if (autoScale){
# # # # #       # Make a basic volcano plot
# # # # #       with(tbl, plot(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, main="Volcano plot",xlab="log Differrential Expression",ylab="-log10(pValue)"))
# # # # #     } else {
# # # # #       # Make a basic volcano plot
# # # # #       with(tbl, plot(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, main="Volcano plot", xlim=xlimit,ylim=ylimit))
# # # # #     }
# # # # # 
# # # # #   # Add colored points: red if padj < 0.05, orange of log2FC > 1, green if both)
# # # # #   with(subset(tbl, p.adjust(Significance_pvalue, "none", length(Significance_pvalue)) <= pValCutoff), points(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, col="red"))
# # # # #   with(subset(tbl, (Value_LogDiffExp <= DiffExpLowerCutOff | Value_LogDiffExp >= DiffExpUpperCutOff)), points(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, col="orange"))
# # # # #   with(subset(tbl, p.adjust(Significance_pvalue, "none", length(Significance_pvalue)) <= pValCutoff & (Value_LogDiffExp <= DiffExpLowerCutOff | Value_LogDiffExp>=DiffExpUpperCutOff)), points(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, col="green"))
# # # # # 
# # # # #   # Label points with the textxy function from the calibrate plot
# # # # #   if (labels){
# # # # # #     library(calibrate)
# # # # #     with(subset(tbl, p.adjust(Significance_pvalue, "none", length(Significance_pvalue)) <= pValCutoff & (Value_LogDiffExp <= DiffExpLowerCutOff | Value_LogDiffExp >= DiffExpUpperCutOff)), textxy(Value_LogDiffExp, -log10(Significance_pvalue), labs=Name_GeneSymbol, cex=cex))
# # # # #   }
# # # # # 
# # # # #   }
# # # # # }

################################################




# mouseSigGSE2988.xls
# path_to_write="/opt/raid10/genomics/mehdi/ilincs/output/"
volcanoPlotly <- function(sigFile=NULL, sigID=NULL, sigType="notENCODE", xlimit=c(-6,6), ylimit=c(0,10), DiffExpLowerCutOff=-2, DiffExpUpperCutOff=2, 
                          pValCutoff=0.05, labels=FALSE, cex=0.6, path_to_write="/mnt/raid/tmp/", output="json", debuging=FALSE)
{
# sigPlot <- function(sigFile=NULL, sigID=NULL, path_to_write="/mnt/raid/tmp/", output="json", debuging=FALSE)
# {
    
  res <- list(fileName=NA, Remark="Done", obj=NA)
  
  if(is.null(sigFile)) {
    if(is.null(sigID)) {
      res$Remark <- "Missing signature!"
      return(res)
    }
  }
  
  pValCutoff <- -log10(pValCutoff)
  
a <- Sys.time()
print(paste("volcanoPlotly called with following parameters at:",a))
print(paste("sigFile:", sigFile))
print(paste("sigID:", sigID))
print(paste("path_to_write:", path_to_write))
print(paste("debuging:", debuging))

if(!is.null(sigFile)) {
  dat <- try(read.table(file=paste0(path_to_write, sigFile), header=TRUE, stringsAsFactors=FALSE,  sep="\t")) # , quote=""
} else {
  dat <- try(downloadSignature(sigID = sigID, path_to_write = path_to_write, display=TRUE, geneNames=TRUE, write=FALSE, debuging=debuging))
}

if(class(dat) != "data.frame") {
    res$Remark <- "No appropriate signature ID, file or path!"
    return(res)
}


if (sigType=="ENCODE") {
  dat <- dat[, c("ID_geneid", "Name_GeneSymbol", "Score", "Probability")]
  res$Remark <- "Wrong sgnature type!"
  return(res)
  
  
} else {
  dat <- dat[, c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")]
  dat[,4] <- -log10(dat[,4])
  
  if(dim(dat)[2] == 0) {
    res$Remark <- "No genes in the signature!"
    return(res)
  }
  
  if(dim(dat)[1] <= 4) {
    res$Remark <- "Not enough columns provided!"
    return(res)
  }
  
  dat$colors <- "black"
  dat$colors[dat$Significance_pvalue >= pValCutoff] <- "red"
  dat$colors[dat$Value_LogDiffExp <= DiffExpLowerCutOff | dat$Value_LogDiffExp >= DiffExpUpperCutOff] <- "orange"
  dat$colors[dat$Significance_pvalue >= pValCutoff & (dat$Value_LogDiffExp <= DiffExpLowerCutOff | dat$Value_LogDiffExp >= DiffExpUpperCutOff)] <- "green"
  colnames(dat)[4] <- "-log10(Significance_pvalue)"
  
  
  if("console" %in% output & "json" %in% output) {
    res$Remark <- "Can not return json and console! ;-)"
    return(res)
  }
  
  if("png" %in% output | "pdf" %in% output | "ggplot" %in% output | "console" %in% output) {
    plt <- ggplot2::ggplot(dat, ggplot2::aes_(x = as.name(names(dat)[3]), y = as.name(names(dat)[4]))) + # , fill=as.name(colnames(dat)[2])) + 
      ggplot2::geom_point(color = dat$colors, size = 0.5) + ggplot2::theme_bw()
    #     res$obj <- plt
  }
  
  
  
  
  sessionID <- generateSID()
  res$fileName <- paste0("volcanoPlot_", sessionID)
  
  if("png" %in% output) {
    file <- paste0(path_to_write, "volcanoPlot_", sessionID, ".png")
    ggplot2::ggsave(file, plot=plt)
    #     try(dev.off())
    # # # 	return(paste0("Done file: ", path_to_write, "gseaPlot_", sessionID, ".png"))
  }
  if("pdf" %in% output) {
    file <- paste0(path_to_write, "volcanoPlot_", sessionID, ".pdf")
    ggplot2::ggsave(file, plot=plt)
    #     try(dev.off())
    # # # 	return(paste0("Done file: ", path_to_write, "gseaPlot_", sessionID, ".pdf"))
  }
  if("ggplot" %in% output) {
    # 	res$fileName <- plt
    save(plt, file = paste0(path_to_write, "volcanoPlot_", sessionID, ".RData"))
    # # # 	return(res)
  }
  
  if("json" %in% output | "widget" %in% output) {
    # # # 	pl <- plotly::ggplotly(plt)
    pal <- c("black", "red", "orange", "green")
    pal <- setNames(pal, pal)
    pl <- plotly::plot_ly(dat, x = ~Value_LogDiffExp, y = ~`-log10(Significance_pvalue)`, text = dat$Name_GeneSymbol, mode = "markers", type = "scatter", 
                          color = ~colors, colors=pal, marker=list(size = 5), showlegend = FALSE)
  }
  
  if("widget" %in% output) {
    htmlwidgets::saveWidget(pl, paste0(path_to_write, "volcanoPlot_", sessionID, ".html"), selfcontained = TRUE)
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
}


return(res)
}

