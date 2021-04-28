#Filter SNP-index variants 
#' @title Filter SNP-index Variants 
#' @description This function filters out variants in the data frame returned by the calc_SNPindex() function 
#' which  do not  fall  between minimum DP, maximum DP, minimum SNP-index, maximum SNP-index, and minimum GQ for both bulks (wild-type and mutant).
#'
#' @param vcf.df.SNPindex vcf dataframe
#' @param min.SNPindex min value allowed for the SNP index (default=0.3)
#' @param max.SNPindex max value allowed for the SNP index (default=0.9)
#' @param min.DP min value allowed for the read depth (default=50)
#' @param max.DP max value allowed for the read depth (default=200)
#' @param min.GQ min value allowed for the genotype quality (default=99)
#'
#' @return data frame containing variant and SNP-index information on variants that passed the filtering process.
#'
#' @details Variants from the data frame returned by the calc_SNPindex() function 
#' are filtered out based on the following criteria:
#' minimum DP, maximum DP, minimum SNP-index, maximum SNP-index, and minimum GQ. 
#'
#' Variants must pass these filters in both bulks to be retained.
#'
#' @importFrom dplyr %>%
#' @export filter_SNPindex
#' @examples
#' ## Default parameters
#' vcf_df_SNPindex_filt <- filter_SNPindex(vcf.df.SNPindex=vcf_df_SNPindex)
#' ## Custom parameters
#' vcf_df_SNPindex_filt <- filter_SNPindex(vcf.df.SNPindex=vcf_df_SNPindex, 
#'                                         min.SNPindex=0.25, 
#'                                         max.SNPindex=0.8, 
#'                                         min.DP=60, 
#'                                         max.DP=250, 
#'                                         min.GQ=98)


filter_SNPindex <- function(vcf.df.SNPindex, min.SNPindex=0.3, max.SNPindex=0.9, min.DP=50, max.DP=200, min.GQ=99){
  
  #Filter by min SNP-index
  vcf.df.SNPindex <- vcf.df.SNPindex %>% dplyr::filter(SNPindex.WT >= min.SNPindex & SNPindex.M >= min.SNPindex &
                                                       SNPindex.WT <= max.SNPindex & SNPindex.M <= max.SNPindex) 
  #Filter by min and max depth
  vcf.df.SNPindex <- vcf.df.SNPindex %>% dplyr::filter(DP.WT >= min.DP & DP.M >= min.DP &
                                                       DP.WT <= max.DP & DP.M <= max.DP) 
  #Filter by min genotype quality
  vcf.df.SNPindex.filt <- vcf.df.SNPindex %>% dplyr::filter(GQ.WT >= min.GQ & GQ.M >= min.GQ) 
  
  return(vcf.df.SNPindex.filt)
}
