#Calculate SNP-ratio

#' @title Calculate SNP-ratio
#' @description This function allows to calculate the SNP-ratio between two bulks (M/WT), as follows:
#'
#' \deqn{(AD_ref.WT/(AD_ref.WT + AD_alt.WT)) - (AD_ref.M/(AD_ref.M + AD_alt.M))}
#'
#' Calculated SNP-ratio values get added in a new column of the data frame (returned by the readBSA_vcf() function).
#' Note that the user can select specific variants to consider, by setting the "Variants" parameter to "SNP" (default) or "all" (respectively referring to SNPs or InDels+SNPs).
#'
#' @references Wachsman et al., 2017
#'
#' @param vcf.df Data frame of the vcf file
#' @param wtBulk Wild-Type pool
#' @param mBulk Mutant pool 
#' @param variants variants to be considered. Default is "SNP" (allowed: "SNP" or "all")
#'
#' @return Data frame containing filtered variant and SNP-ratio information.
#'
#' @details The data frame returned by readBSA_vcf() is filtered by bulk ID to create two separate data frames: 
#' one specific to the wild-type bulk variants information and other one specific to the mutant bulk variants information.
#' In the case the type of variants is set to "SNP", variants corresponding to InDels are discarded by removing the rows in the data frame containing more than three characters in the GT_alleles column (e.g.,"T/AT"corresponds to an insertion and contains 4 characters). 
#' This also applies for deletions, where one of the three characters is a "*". Then, both data frames (wild-type bulk and mutant bulk) are joined by rows having the same value in the chromosome and position columns. 
#' The SNP-ratio is then calculated for each of the variants and the values are added in a new column of the joint data frame, which will be returned by the function.
#' 
#' @importFrom dplyr %>%
#' @export calc_SNPratio
#' @examples
#' ## Calculate SNP-index for both bulks (only SNPs will be considered)
#' vcf_df_SNPratio <- calc_SNPratio(vcf.df=vcf_list$df, 
#'                                  wtBulk="pool_S3781_minus", 
#'                                  mBulk="pool_S3781_plus") 
#' ## Calculate SNP-index considering both InDels and SNPs
#' vcf_df_SNPratio <- calc_SNPratio(vcf.df=vcf_list$df, 
#'                                  wtBulk="pool_S3781_minus", 
#'                                  mBulk="pool_S3781_plus", 
#'                                  variants="all")


calc_SNPratio <- function(vcf.df, wtBulk, mBulk, variants="SNP") {
  
  #Create data frame for each bulk 
  vcf.df.WTbulk <- vcf.df %>% dplyr::filter(Indiv==wtBulk) 
  vcf.df.Mbulk <- vcf.df %>% dplyr::filter(Indiv==mBulk) 
  
  #Stop the program and show message if user selects a not allowed 'variants' value
  if (variants!="SNP" & variants!="all") {
    stop("The allowed values for the 'variants' argument are: 'SNP' or 'all'. The latter will consider both InDels and SNPs.")
  }
  
  #Remove from the data frames rows corresponding to InDel variants (if variants=="SNP")
  if (variants=="SNP") {
    vcf.df.WTbulk <- vcf.df.WTbulk %>% dplyr::filter(nchar(GT_alleles)==3 &
                                                     grepl("[^*]{3}", GT_alleles))
    vcf.df.Mbulk <- vcf.df.Mbulk %>% dplyr::filter(nchar(GT_alleles)==3 &
                                                   grepl("[^*]{3}", GT_alleles))
  }
  
  #Join data frames of each bulk by same chromosome (ChromKey) and position (POS)
  vcf.df.jointBulks <- dplyr::inner_join(vcf.df.WTbulk, 
                                         vcf.df.Mbulk, 
                                         by = c("ChromKey","POS"), 
                                         copy = F, 
                                         suffix = c(".WT", ".M"))
  
  #Calculate SNP-ratio between the two bulks
  vcf.df.SNPratio <- vcf.df.jointBulks %>% dplyr::mutate("SNPratio" = ((as.numeric(AD_ref.WT)/(as.numeric(AD_ref.WT) + as.numeric(AD_alt.WT))) - (as.numeric(AD_ref.M)/(as.numeric(AD_ref.M) + as.numeric(AD_alt.M)))))
  
  return(vcf.df.SNPratio)
}
