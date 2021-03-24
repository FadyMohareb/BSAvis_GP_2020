#Filter variants 
#'
#' This function allows to filter variants.
#'
#' @param vcf.df.bulk
#' @param min_SNPindex min value allowed for the SNP index (default=0.3)
#' @param max_SNPindex max value allowed for the SNP index (default=0.9)
#' @param min_DP min value allowed for the read depth (default=50)
#' @param max_DP max value allowed for the read depth (default=200)
#' @param min_GQ min value allowed for the genotype quality (default=99)
#'
#' @return Filtered bulks (type: data frame)
#'
#' @export


filter_variants <- function(vcf.df.bulks, min_SNPindex=0.3, max_SNPindex=0.9, min_DP=50, max_DP=200, min_GQ=99){
  
  #Filter by min SNP-index
  vcf.df.bulks <- vcf.df.bulks %>% dplyr::filter(SNPindex.WT >= min_SNPindex & SNPindex.M >= min_SNPindex &
                                                 SNPindex.WT <= max_SNPindex & SNPindex.M <= max_SNPindex) 
  #Filter by min and max depth
  vcf.df.bulks <- vcf.df.bulks %>% dplyr::filter(DP.WT >= min_DP & DP.M >= min_DP &
                                                 DP.WT <= max_DP & DP.M <= max_DP) 
  #Filter by min genotype quality
  vcf.df.bulks.filt <- vcf.df.bulks %>% dplyr::filter(GQ.WT >= min_GQ & GQ.M >= min_GQ) 
  
  return(vcf.df.bulks.filt)
}