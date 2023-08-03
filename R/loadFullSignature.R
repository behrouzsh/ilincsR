
#' A function to retrieve a previously created full signature 
#'
#' This function allows you to load/retrieve a full signature which is created in the initial step of analysis.
#' @param sessionID the initial session ID.
#' @param path_to_write This optional parameter specifies the path to load the signature from. It should be the same as the path which signature is saved in.
#' @param Value The output of this function is in the form of a data.frame. First column is geneIDs, second column is geneNames, 
#'	third column is coefficients and the last column is their probability values for each gene from statistical analysis (ttest).
#' @param authors M. Fazel-Najafabadi
#' @keywords SignatureTable.
#' @export 
#' @examples
#' ## not run
#' fullSig <- loadFullSignature("Fri_Jun_9_15_21_03_2017_4119651")
#' ## end not run

loadFullSignature <- function(sessionID, path_to_write="/mnt/raid/tmp/") {
  if (file.exists(paste(path_to_write, "temp_volcano_", sessionID, ".RData", sep=""))) {
      print(paste("loading . . . "))
      load(paste(path_to_write, "temp_volcano_", sessionID, ".RData", sep=""))
      
      SignatureTable <- SignatureTable[, c("PROBE", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")]
      return(SignatureTable)
  } else { 
      print(paste("File not found!"))
      return("Full signature not found.")
  }
}
