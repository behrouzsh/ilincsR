
#' A function to create the concordance table between a signature and a library.
#'
#' This function allows you to ...
#' @param data1 A numeric vector as the input signature.
#' @param data2 A numeric vector or a matrix in the same row orther as in data1.
#' @param cortype The correlation type like "pearson" or "spearman".
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @export 
#' @examples
#' res <- cor_with_pval(data1, data2)

cor_with_pval <- function(data1, data2, cortype="pearson", alternative="two.sided", topSigs=FALSE) {
## data1 and data2 should be data.frames or named vectors. 

    if(!is.matrix(data1)) data1 <- as.matrix(data1)
    if(!is.matrix(data2)) data2 <- as.matrix(data2)
    if (!is.numeric(data1) | !is.numeric(data2)) stop("data1 and data2 should be numeric (vector or matrix)")
    if (dim(data1)[2] > 1 & dim(data2)[2] > 1) stop("One of data1 or data2 should be a vector not matrix)")
    if (dim(data1)[2] == 1) { colnames(data1)[1] <- "data1" }
    if (dim(data2)[2] == 1) { colnames(data2)[1] <- "data2" }
    if (dim(data1)[1] <= 2 | dim(data2)[1] <= 2) return(data.frame(signatureID=NA, similarity=NA, PValue=NA, nGenes=NA))
# # # # # d1 <- which(apply(data1, 2, function(i) length(i[!is.na(i)])) > 2)
# # # # # data1 <- data1[, d1, drop=FALSE]
# # # # # d2 <- which(apply(data2, 2, function(i) length(i[!is.na(i)])) > 2)
# # # # # data2 <- data2[, d2, drop=FALSE]

    if (dim(data1)[2] > 1 & dim(data2)[2] == 1){
	    dat <- data1
	    data1 <- data2
	    data2 <- dat
    }

# # #     tmp <- as.data.frame(data.table::rbindlist(tmp))
# # #     tmp <- cbind(colnames(y), tmp, stringsAsFactors=FALSE)
# # # #     n <- t(!is.na(x)) %*% (!is.na(y))
# # #     colnames(tmp) <- c("signatureID", "Similarity", "PValue", "nGenes")
# # #     tmp[,4] <- tmp[,4] + 2
# # #     tmp
# # # }

    tmp <- cor(data1, data2, use="pairwise.complete.obs", method=cortype)
    df <- (t(!is.na(data1)) %*% (!is.na(data2))) - 2
    df[df < 0] <- NaN
    res <- data.frame(signatureID=colnames(tmp), Similarity=tmp[1,], PValue=NA, nGenes=as.integer(df[1,]), stringsAsFactors=FALSE)
    if(topSigs) {
	bp <- boxplot(res$Similarity, plot = FALSE)
	res <- res[res[,2] <= bp$stats[1,1] | res[,2] >= bp$stats[5,1],]
    }
    
#     stat <- sqrt(df[1,]) * tmp[1,] / sqrt(1 - tmp[1,] * tmp[1,])
#     tp <- rbind(pt(stat, df[1,]), pt(stat, df[1,], lower.tail=FALSE))
    stat <- sqrt(res[,4]) * res[,2] / sqrt(1 - res[,2] * res[,2])
    tp <- rbind(pt(stat, res[,4]), pt(stat, res[,4], lower.tail=FALSE))
    
    pv <- switch(alternative,
		       "less" = tp[1,],
		       "greater" = tp[2,],
		       "two.sided" = 2 * apply(tp, 2, min)
		)

    res[,4] <- res[,4] + 2
    res[,3] <- pv
#     res[,3] <- as.numeric(formatC(res[,3], format = "e", digits = 4))
    
return(res)
}

