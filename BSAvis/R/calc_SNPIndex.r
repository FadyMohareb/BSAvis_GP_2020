#Calculate SNP-index

#' @title Calculate SNP-index
#' @description This function calculates, for both wild-type and mutant bulks, the SNP-index value of each variant stored inside the data frame (returned by the readBSA_vcf() function), 
#' applying the following formula:
#' \deqn{SNPindex=AD_alt/(AD_ref + AD_alt)} 
#' Bulks get joined together in a single data frame.
#' 
#' Note that the user can select specific variants to consider, 
#' by setting the "variants" parameter to "SNP" (default) or "all" (InDels+SNPs).
#'
#' @param vcf.df Data frame of the vcf file
#' @param wtBulk Wild-Type pool
#' @param mBulk Mutant pool 
#' @param variants variants to be considered. Default is "SNP" (allowed: "SNP" or "all")
#'
#' @return Data frame containing variand and SNP-index information for each bulk.
#'
#' @details The data frame returned by readBSA_vcf() is filtered by bulk ID to create two separate data frames: 
#' one specific to the wild-type bulk variants information and other one specific to the mutant bulk variants information.
#' For every variant in each data frame, the SNP-index value gets calculated and added as a new column.
#' 
#' In the instance that the type of variants to be considered is set to "SNP",the variants corresponding to InDel sare discarded.  
#' The rows containing more than three characters in the GT_alleles column (e.g.,"T/AT"corresponds to an insertion and contains 4 characters) 
#' or containing three characters but one of them being "*", meaning a deletion, are removed from the data frame. 
#' 
#' Both dataframes are then joined by same chromosome and position to ensure that all the merged rows contain information on the same genomic position.
#' This data frame is returned by the function
#'
#' @export calc_SNPindex
#' @examples
#' ## Calculate SNP-index for both bulks (only SNPs will be considered)
#' vcf_df_SNPindex <- calc_SNPindex(vcf.df=vcf_list$df, 
#'                                  wtBulk="pool_S3781_minus", 
#'                                  mBulk="pool_S3781_plus") 
#'
#' ## Calculate SNP-index considering both InDels and SNPs
#' vcf_df_SNPindex <- calc_SNPindex(vcf.df=vcf_list$df, 
#'                                  wtBulk="pool_S3781_minus", 
#'                                  mBulk="pool_S3781_plus",
#'                                  variants="all") 


calc_SNPindex <- function(vcf.df, wtBulk, mBulk, variants="SNP") {
  
  #Create data frame for each bulk AND include a new column with SNP index
  vcf.df.wtBulk <- vcf.df %>% dplyr::filter(Indiv==wtBulk) %>% dplyr::mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))
  vcf.df.mBulk <- vcf.df %>% dplyr::filter(Indiv==mBulk) %>% dplyr::mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))
  
  #Stop the program and show message if user selects a non-allowed 'variants' value
  if (variants!="SNP" & variants!="all") {
    stop("The allowed values for the 'variants' argument are: 'SNP' or 'all'. The latter will consider both InDels and SNPs.")
  }
  
  #Remove from data frames those rows corresponding to InDel variants (if variants=="SNP")
  if (variants=="SNP") {
    vcf.df.wtBulk <- vcf.df.wtBulk %>% dplyr::filter(nchar(GT_alleles)==3 &
                                                     grepl("[^*]{3}", GT_alleles))
    vcf.df.mBulk <- vcf.df.mBulk %>% dplyr::filter(nchar(GT_alleles)==3 &
                                                   grepl("[^*]{3}", GT_alleles))
  }
  
  #Join data frames of each bulk by same chromosome (ChromKey) and position
  vcf.df.SNPindex <- dplyr::inner_join(vcf.df.wtBulk, 
                                       vcf.df.mBulk, 
                                       by = c("ChromKey","POS"), 
                                       copy = F, 
                                       suffix = c(".WT", ".M"))
  
  return(vcf.df.SNPindex)
}

