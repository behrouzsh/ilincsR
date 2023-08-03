
#' A function for ...
#'
#' This function allows you to ...
#' @param sigfile ...
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



volcanoPlot <- function(sigfile, sigType="notENCODE", xlimit=c(-6,6), ylimit=c(0,10), DiffExpLowerCutOff=-2, DiffExpUpperCutOff=2, 
			pValCutoff=0.05, labels=FALSE, cex=0.6, path_to_write="/mnt/raid/tmp/")
{

textxy <- function (X, Y, labs, m = c(0, 0), cex = 0.5, offset = 0.8, ...)
{
    posXposY <- ((X >= m[1]) & ((Y >= m[2])))
    posXnegY <- ((X >= m[1]) & ((Y < m[2])))
    negXposY <- ((X < m[1]) & ((Y >= m[2])))
    negXnegY <- ((X < m[1]) & ((Y < m[2])))
    if (sum(posXposY) > 0)
        text(X[posXposY], Y[posXposY], labs[posXposY], adj = c(0.5 - offset, 0.5 - offset), cex = cex, ...)
    if (sum(posXnegY) > 0)
        text(X[posXnegY], Y[posXnegY], labs[posXnegY], adj = c(0.5 - offset, 0.5 + offset), cex = cex, ...)
    if (sum(negXposY) > 0)
        text(X[negXposY], Y[negXposY], labs[negXposY], adj = c(0.5 + offset, 0.5 - offset), cex = cex, ...)
    if (sum(negXnegY) > 0)
        text(X[negXnegY], Y[negXnegY], labs[negXnegY], adj = c(0.5 + offset, 0.5 + offset), cex = cex, ...)
}

  tbl <- read.table(paste0(path_to_write, sigfile), header=TRUE, sep="\t", stringsAsFactors=FALSE) # , quote = ""
  
  autoScale <- TRUE
  if (sigType=="ENCODE") {

    if (autoScale){
      # Make a basic volcano plot
      with(tbl, plot(Score, Probability, pch=20, main="Scatter plot")) #ENC
    } else {
      # Make a basic volcano plot  
      with(tbl, plot(Score, Probability, pch=20, main="Volcano plot", xlim=xlimit,ylim=ylimit)) #ENC      
    }


    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(tbl, p.adjust(Probability, "none", length(Probability)) < pValCutoff), points(Score, Probability, pch=20, col="red"))
    with(subset(tbl, Score >= DiffExpUpperCutOff), points(Score, Probability, pch=20, col="orange"))
    with(subset(tbl, p.adjust(Probability,"none", length(Probability)) <= pValCutoff & (Score >= DiffExpUpperCutOff)), points(Score, Probability, pch=20, col="green"))
    
    # Label points with the textxy function from the calibrate plot
    if (labels){
#       library(calibrate)
      with(subset(tbl, p.adjust(Probability,"none",length(Probability)) <= pValCutoff & (Score >= DiffExpUpperCutOff)), textxy(Score, Probability, labs=Name_GeneSymbol, cex=cex))
    }

  } else {

#Example
    if (autoScale){
      # Make a basic volcano plot
      with(tbl, plot(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, main="Volcano plot",xlab="log Differrential Expression",ylab="-log10(pValue)"))
    } else {
      # Make a basic volcano plot
      with(tbl, plot(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, main="Volcano plot", xlim=xlimit,ylim=ylimit))
    }

  # Add colored points: red if padj < 0.05, orange of log2FC > 1, green if both)
  with(subset(tbl, p.adjust(Significance_pvalue, "none", length(Significance_pvalue)) <= pValCutoff), points(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, col="red"))
  with(subset(tbl, (Value_LogDiffExp <= DiffExpLowerCutOff | Value_LogDiffExp >= DiffExpUpperCutOff)), points(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, col="orange"))
  with(subset(tbl, p.adjust(Significance_pvalue, "none", length(Significance_pvalue)) <= pValCutoff & (Value_LogDiffExp <= DiffExpLowerCutOff | Value_LogDiffExp>=DiffExpUpperCutOff)), points(Value_LogDiffExp, -log10(Significance_pvalue), pch=20, col="green"))

  # Label points with the textxy function from the calibrate plot
  if (labels){
#     library(calibrate)
    with(subset(tbl, p.adjust(Significance_pvalue, "none", length(Significance_pvalue)) <= pValCutoff & (Value_LogDiffExp <= DiffExpLowerCutOff | Value_LogDiffExp >= DiffExpUpperCutOff)), textxy(Value_LogDiffExp, -log10(Significance_pvalue), labs=Name_GeneSymbol, cex=cex))
  }

  }
}
