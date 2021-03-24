#Calculate delta(SNP Index)
#'
#' This function is used to calculate the delta(SNP index) using the calculated mean of the M and WT indices, as follows:
#' \deqn{mean_SNPindex.M - mean_SNPindex.WT)}. 
#'
#' @param vcf.df.window data frame obtained after running the slidingWindow() function
#'
#' @return data frame containing the delta(SNP index) values
#'
#' @export
#' @examples
#' calc_deltaSNPindex()


calc_deltaSNPindex <- function(vcf.df.window) {
  
  #Calculate delta SNP-index and add it as new column
  vcf.df.deltaSNPindex <- vcf.df.window %>% dplyr::mutate("delta_SNPindex" = (mean_SNPindex.M - mean_SNPindex.WT))
  return(vcf.df.deltaSNPindex)
}