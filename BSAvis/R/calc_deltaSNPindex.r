#Calculate delta(SNP Index)

#' @title Calculate delta(SNP index)
#' @description This function is used to calculate the ∆(SNP-index) by subtracting the SNP-index value of the wild-type bulk from the SNP-index value of the mutant bulk. 
#'
#' \deqn{mean_M_SNPindex - mean_WT_SNPindex.WT}
#'
#' @param SNPindex.windows data frame obtained after running the slidingWindow()function
#'
#' @return Data frame containing start, mid and stop positions of each window in the chromosome, as well as the corresponding mean SNP-index value for each bulk and ∆(SNP-index)value.
#'
#' @details For each row in the data frame returned by slidingWindow(), corresponding to a different window, ∆(SNP-index) is calculated by subtracting the SNP-index value of the reference bulk to the SNP-index value of the mutant bulk. 
#' The values of ∆(SNP-index)ineach window are added to the data frame as a new column.
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
