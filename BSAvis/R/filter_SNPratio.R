#Filter SNP-ratio variants 
#' @title Filter SNP-index Variants 
#' @description This functions allows to filter out bulk variants stored inside the data frame (returned by the calc_SNPratio() function) which do not fall between the 
#' minimumSNP-index, minimum DP and maximum DP. 
#'
#' @param vcf.df.SNPratio vcf dataframe
#' @param min.SNPratio min value allowed for the SNP index (default=0.1)
#' @param min.DP min value allowed for the read depth (default=50)
#' @param max.DP max value allowed for the read depth (default=200)
#'
#' @return Data frame containing filtered variant information.
#'
#' @details Variants in the data frame returned by calc_SNPratio() with SNP-ratio values less below the minimum SNP-ratio value, as well as variants which do not fall between the given (or default) minimum and maximum DP values in both bulks, are discarded and removed from the final data frame returned by the function.
#'
#' @importFrom dplyr %>%
#' @export filter_SNPratio
#' @examples
#' ## Default parameters
#' vcf_df_SNPratio_filt <- filter_SNPratio(vcf.df.SNPratio=vcf_df_SNPratio)
#' ## Custom parameters
#' vcf_df_SNPratio_filt <- filter_SNPratio(vcf.df.SNPratio=vcf_df_SNPratio, 
#'                                         min.SNPratio=0.3, 
#'                                         min.DP=60, 
#'                                         max.DP=250)


filter_SNPratio <- function(vcf.df.SNPratio, min.SNPratio=0.1, min.DP=50, max.DP=200){
  #Filter by min SNP-ratio
  vcf.df.SNPratio <- vcf.df.SNPratio %>% dplyr::filter(SNPratio >= min.SNPratio) 
  
  #Filter by min and max depth
  vcf.df.SNPratio.filt <- vcf.df.SNPratio %>% dplyr::filter(DP.WT >= min.DP & DP.M >= min.DP &
                                                            DP.WT <= max.DP & DP.M <= max.DP) 
  
  return(vcf.df.SNPratio.filt)
}