
#' A function for ...
#'
#' This function allows you to ...
#' @param sigfile This should be a signature file created or uploaded by user.
#' @param sigID This is a precomputed signature from any sig-library.
#' @param glist A list of Entrez gene IDs to enrich.
#' @param path_to_write This parameter specifies the path which user wants to save the results.
#' @param output Type of the output plot including "png", "pdf", "ggplot", "json" and "widget".
#' @param debuging This is for database testing purposes only. The default database server used in all work flows is "eh3" but it can be set to 
#'	alternative database server "gimm2" if debuging=TRUE dev servers will be used.
#' @param authors M. Fazel-Najafabadi
#' @keywords filter
#' @export
#' @examples
#' ## Do not run
#' sigfile <- "sig_Thu_Apr_12_12_45_57_2018_298226.xls"
#' sigID <- "LINCSCP_109"
#' path_to_write <- "/opt/raid10/genomics/mehdi/ilincs/output/"
#' plt <- sigPlot(sigfile, sigID, path_to_write, debuging=TRUE, output=TRUE)
#' ## End do not run


gseaPlot <- function(sigFile=NULL, sigID=NULL, glist, path_to_write="/mnt/raid/tmp/", output=c("png", "pdf"), fontSize=15, debuging=FALSE)
{


# # source("https://bioconductor.org/biocLite.R")
# library(BiocInstaller)
# # biocLite("fgsea")
# library(fgsea)
# data(examplePathways)
# data(exampleRanks)
# fgseaRes <- fgsea(pathways = examplePathways, stats = exampleRanks, minSize=15, maxSize=500, nperm=100000)
# head(fgseaRes[order(pval), ])
# sum(fgseaRes[, padj < 0.01])
# 
# topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# 
# pdf("/mnt/raid/tmp/fgsea2.pdf")
# plotEnrichment(examplePathways[[head(fgseaRes[order(pval), ], 1)$pathway]], exampleRanks)# + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)
# plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, gseaParam = 0.5)
# dev.off()
# 
# sig <- downloadSignature("LINCSCP_1254",write=F,debuging = T)
# sts <- sig$signature
# sts <- sts[order(sts$Value_LogDiffExp),]
# stats <- sts$Value_LogDiffExp
# names(stats) <- sts$ID_geneid
# glist <- sample(sts$ID_geneid,120)
# plotEnrichment(glist, stats)
# pl <- ggplotly(plotEnrichment(glist, stats))

###################
{
function (pathway, stats, gseaParam = 1, ticksSize = 0.2, fontSize=13) 
{

  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  sp <- stats[ord]
  sp <- rev(split(sp,cut(sp, seq(-20, 20, 0.2), include.lowest=TRUE, labels=FALSE)))
  xx <- data.frame(Expr_Data=sapply(sp, median), Bin=sapply(sp, length))
  xx$high <- cumsum(xx$Bin)
  xx$low <- xx$high - xx$Bin
  g <- ggplot2::ggplot(toPlot, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color = "green", size = 0.1) + 
    ggplot2::geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    ggplot2::geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
    ggplot2::geom_hline(yintercept = 0, colour = "black") + 
    ggplot2::geom_line(color = "green") + 
    ggplot2::theme_bw() + 
    ggplot2::geom_segment(data = data.frame(x = pathway), mapping = ggplot2::aes(x = x, y = diff/2, xend = x, yend = 0), size = ticksSize) + 
    ggplot2::geom_rect(data=xx, ggplot2::aes(x=NULL, y=NULL, xmin = low, xmax = high, ymin = min(toPlot$y)-diff/2, ymax = min(toPlot$y), fill = Expr_Data)) + 
    ggplot2::scale_fill_gradient2(low = "blue",mid = "white",high = "red") + 
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) + 
    ggplot2::theme(text = ggplot2::element_text(size = fontSize), axis.text = ggplot2::element_text(size = fontSize)) +
    ggplot2::labs(x = "rank", y = "enrichment score")
  g
}
} -> pltEnrichment

# sig <- downloadSignature("LINCSCP_1254",write=F,debuging = T)
# sts <- sig$signature
# sts <- sts[order(sts$Value_LogDiffExp),]
# stats <- sts$Value_LogDiffExp
# names(stats) <- sts$ID_geneid
# glist <- sample(sts$ID_geneid,120)
# plotEnrichment(glist, stats)
# pl <- ggplotly(plotEnrichment(glist, stats))
    res <- list(fileName=NA, Remark="Done", obj=NA)
    if("console" %in% output & "json" %in% output) {
	res$Remark <- "Can not return json and console! ;-)"
	return(res)
    }
    if(!is.null(sigFile)) {
	sig <- uploadedSigProcess(exp=sigFile, userUploadedProfile=path_to_write, write=FALSE, debuging=debuging)$Sig
	sig <- sig[,-wich(colnames(sig)=="PROBE")]
	colnames(sig) <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
    } else if(!is.null(sigID)) {
	sig <- downloadSignature(sigID, write = FALSE, debuging = debuging)
	sig <- sig[,-1]
	sig <- sig[,-which(colnames(sig)=="PROBE")]
    } else {
	res$Remark <- "Missing signature!"
	return(res)
    }
    sig <- sig[!is.na(sig[,3]),]
    sig <- sig[order(sig[,3], decreasing=TRUE),]
    stats <- sig[,3]
    names(stats) <- sig$Name_GeneSymbol # ID_geneid
    if(sum(is.na(as.integer(glist))) > length(glist)/2) {
	glist <- sig$Name_GeneSymbol[na.omit(match(glist, sig$Name_GeneSymbol))] # na.omit(match(glist, names(stats)))
    } else {
	glist <- sig$Name_GeneSymbol[na.omit(match(glist, sig$ID_geneid))] # na.omit(match(glist, names(stats)))
    }
    plt <- pltEnrichment(glist, stats, fontSize=fontSize)

    sessionID <- generateSID()
    res$fileName <- paste0("gseaPlot_", sessionID)
    if("png" %in% output) {
	file <- paste0(path_to_write, "gseaPlot_", sessionID, ".png")
	ggplot2::ggsave(file, plot=plt, dpi=600)
    }
    if("pdf" %in% output) {
	file <- paste0(path_to_write, "gseaPlot_", sessionID, ".pdf")
	ggplot2::ggsave(file, plot=plt, dpi=600)
    }
    if("ggplot" %in% output) {
	save(plt, file = paste0(path_to_write, "gseaPlot_", sessionID, ".RData"))
    }
    if("json" %in% output | "widget" %in% output) {
	pl <- plotly::ggplotly(plt)
	x1 <- pl$x$data[[1]]$x
	x1[2:(length(x1)-1)] <- names(stats)[x1[2:(length(x1)-1)]]
	x5 <- pl$x$data[[5]]$x
	x5 <- names(stats)[x5]
	pl$x$data[[1]]$text <- paste0("Gene Name: ", x1, "<br />Enrichment Score: ",pl$x$data[[1]]$y)
	pl$x$data[[5]]$text <- paste0("Gene Name: ", x5) # , "<br />y: ",pl$x$data[[5]]$y)
    }
    if("widget" %in% output) {
	htmlwidgets::saveWidget(pl, paste0(path_to_write, "gseaPlot_", sessionID, ".html"), selfcontained = TRUE)
    }
    if("console" %in% output) {
	res$obj <- plt
    }
    if("json" %in% output) {
	res$obj <- plotly::plotly_json(pl, FALSE)
    }
    if(sum(output %in% c("png", "pdf", "widget", "json", "ggplot")) == 0) {
	res$Remark <- "Error!!!"
# 	return(res)
    }
    return(res)
	
}
