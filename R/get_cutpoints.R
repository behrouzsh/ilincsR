
#' A function to return percentiles of different signature library concordances
#'
#' This function allows you to get the 0.5, 1.0, 99.0 and 99.5 percentiles of different concordance signature libraries.
#' @param libName The name of signature library user wants to get percentiles for. The default is all.
#' @param level any of 0.5, 1.0, 99.0 and 99.5 percentiles equivalent to positive 0.5 and 1.0 percentiles and negative 0.5 and 1.0 percentiles.
#`		Values are c("neg0.5", "neg1.0", "pos1.0", "pos0.5") The default is all.
#' @param author M. Fazel-Najafabadi
#' @keywords Concordance.
#' @export 
#' @examples
#' ## not run
#' get_cutpoints()
#' # returns all available percentiles for all library concordances.
#' ## end not run

get_cutpoints <- function(libName=c("LIB_1", "LIB_2", "LIB_3", "LIB_5", "LIB_6", "LIB_7", "LIB_8", "LIB_9", "LIB_10", "LIB_11", "LIB_12"), 
			  level=rev(c("pos05", "pos10", "neg10", "neg05"))) 
{

    cutoffs <- data.frame(
		    libName=c("LIB_1", 		"LIB_2", 	"LIB_3", 	"LIB_5", 	"LIB_6", 	"LIB_8", 	"LIB_9", 	"LIB_10", 	"LIB_11", 	"LIB_12", 	"LIB_13", 	"LIB_14"),
		    "pos05"=c(0.176644, 	0.259328, 	8.946850, 	0.321475, 	0.364163, 	0.103731, 	0.110531, 	0.152343, 	0.223605, 	0.112591, 	0.124612, 	0.105921),
		    "pos10"=c(0.15035, 		0.212926, 	4.669846, 	0.28295, 	0.322974, 	0.090643, 	0.060243, 	0.107545, 	0.225157, 	0.101259, 	0.105130, 	0.111258),
		    "neg10"=c(-0.145081, 	-0.125341, 	-7.238341, 	-0.199421, 	-0.224126, 	-0.119243, 	-0.105443, 	-0.104525, 	-0.234185, 	-0.063581, 	-0.184801, 	-0.055801),
		    "neg05"=c(-0.16911, 	-0.144081, 	-11.648286, 	-0.234868, 	-0.262584, 	-0.136783, 	-0.120483, 	-0.115977, 	-0.236064, 	-0.115476, 	-0.130063, 	-0.109466),
		    stringsAsFactors=FALSE
		    )
res <- cutoffs[cutoffs$libName %in% libName, c("libName", level)]
rownames(res) <- NULL
return(res)
}

