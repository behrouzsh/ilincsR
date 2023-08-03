
#' A function to calculate weighted concordance between a signature and a library.
#'
#' This function allows you to ...
#' @param data1 A numeric vector as the input signature.
#' @param data2 A numeric vector or a matrix in the same row orther as in data1.
#' @param weight A numeric vector or a matrix in the same row orther as in data1 
#' and data2 and the same number of columns in data2.
#' @param authors M. Fazel-Najafabadi
#' @keywords ...
#' @useDynLib ilincsR mehdi_wt_cor
#' @export 
#' @examples
#' res <- mehdi_wtcor(data1, data2, weight)


mehdi_wtcor <- function(data1, data2=NULL, weight=NULL, alternative="two.sided", pvalue=TRUE, topSigs=FALSE){
  q <- as.matrix(data1)
  if(ncol(q) > 1) stop("first argument should be a vector or one column array!")
  if(is.null(data2)){
    data2 <- data1
  }
  if(is.null(weight)){
    warning("Missing weights, using equal weights for all variables!")
    weight <- matrix(rep(1, length(data2)), nrow=dim(q)[1])
  }
  if(any(weight < 0, na.rm=TRUE)) stop("Weights should be all positive!")
  if(any(weight == Inf)) stop("Weights contain Inf, not allowed!")
  if(length(weight[is.na(weight)]) > 0) {
      warning("Weights contain NA, using (0) instead!")
      weight[is.na(weight)] <- 0
  }
  if(is.matrix(data2)) r <- data2 else r <- as.matrix(data2)
  if(is.matrix(weight)) w <- weight else w <- as.matrix(weight)
  if(ncol(r) != ncol(w)) stop("Only one column weight is allowed for each pair of data!")
  if((nrow(q) != nrow(w)) | (nrow(r) != nrow(w))) stop("All data and weights should have same number and order of rows!")
  out <- .Call("mehdi_wt_cor", q, r, w, NAOK=TRUE, PACKAGE="ilincsR")
  out[2,out[1,]>=1] <- Inf
  out[2,out[1,]<=-1] <- -Inf
  ## C code for this function is provided by Mehdi Fazel-Najafabadi
  if(!is.null(colnames(r))) colnames(out) <- colnames(r)
  out <- out[,!is.na(out[1,]),drop=FALSE]
  
  if(topSigs) {
      bp <- boxplot(out[1,], plot = FALSE)
      tp <- out[1,] <= bp$stats[1,1] | out[1,] >= bp$stats[5,1]
      out <- out[,tp,drop=FALSE]
  }
  if(pvalue) {
      tp <- rbind(pt(out[2,], out[3,]), pt(out[2,], out[3,], lower.tail=FALSE))
      pv <- switch(alternative,
			  "less" = tp[1,],
			  "greater" = tp[2,],
			  "two.sided" = 2 * apply(tp, 2, min)
		    )
      out <- rbind(out, pv)
      rownames(out) <- c("estimate", "statistic", "df", "pvalue")
  } else {
      rownames(out) <- c("estimate", "statistic", "df")
  }
  
  out
}

