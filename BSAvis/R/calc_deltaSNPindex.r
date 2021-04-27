#Calculate delta(SNP Index)

#' @title Calculate delta(SNP index)
#' @description This function is used to calculate the delta(SNP index) using the calculated mean of the mutant (M) and wild-type (WT) indices, as follows:
#' \deqn{mean_M_SNPindex - mean_WT_SNPindex.WT}
#'
#' @param SNPindex.windows data frame obtained after running the slidingWindow()function
#'
#' @return data frame containing the delta(SNP-index) means values
#'
#' @importFrom dplyr %>%
#' @export calc_deltaSNPindex
#' @examples
#' deltaSNPindex_windows <- calc_deltaSNPindex(SNPindex.windows=SNPindex_windows) 
 

calc_deltaSNPindex <- function(SNPindex.windows) {
  
  #Calculate delta SNP-index and add it as new column
  deltaSNPindex.windows <- SNPindex.windows %>% dplyr::mutate("delta_SNPindex" = (mean_SNPindex.M - mean_SNPindex.WT))
  return(deltaSNPindex.windows)
}
