#Filter SNP-ratio variants 
#' @title Filter SNP-index Variants 
#' @description This function allows to filter variants for the SNP-ratio method.
#'
#' @param vcf.df.SNPratio vcf dataframe
#' @param min.SNPratio min value allowed for the SNP index (default=0.1)
#' @param min.DP min value allowed for the read depth (default=50)
#' @param max.DP max value allowed for the read depth (default=200)
#'
#' @return Filtered bulks for the SNP-ratio method (type: data frame)
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