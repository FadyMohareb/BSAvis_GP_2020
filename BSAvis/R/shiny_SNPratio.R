#Shiny: SNP-ratio wrapper function 

#' @title Shiny SNP-ratio Wrapper Function
#' @description Please note that this function is not to be run manually.
#' Wrapper function for SNP-ratio plots specific to the R-shiny application. 
#' Calls all functions involved in plotting the SNP-ratio values across chromosome positions.
#'
#' @param vcf.list object containing meta information and vcf data frame
#' @param wtBulk Wild-Type pool
#' @param mBulk Mutant pool 
#' @param variants variants to be considered. Default is "SNP" (allowed: "SNP" or "all")
#' @param min.SNPratio min value allowed for the SNP index (default=0.3)
#' @param min.DP min value allowed for the read depth (default=50)
#' @param max.DP max value allowed for the read depth (default=200)
#' @param chrID chromosome ID of interest
#' @param chr chromosome name printed on the plot
#' @param degree LOESS smoothing degree (default=2)
#' @param span LOESS smoothing span (default=0.07)
#' @param ranges axes ranges (x,y)
#' 
#' @details Wrapper function for generting SNP-ratio plots inside the R-Shiny application.
#' The involved functions are the following (in oder): 
#' calc_SNPratio(), filter_SNPratio(), extract_chrIDs() and shinyPlot_SNPratio().
#' 
#' This function is not intended to be run by the user, as it was specifically created to be used by the BSAvis_shiny() function.
#' @importFrom dplyr %>%
#' @export shiny_SNPratio



shiny_SNPratio <- function(vcf.list, wtBulk, mBulk, variants,
                           min.SNPratio, min.DP, max.DP,
                           chrID, chr, degree, span, ranges){
  
  #Calculate SNP-ratio of each variant
  vcf_df_SNPratio <- BSAvis::calc_SNPratio(vcf.list$df, wtBulk, mBulk, variants)
  
  #Filter variants
  vcf_df_SNPratio_filt <- BSAvis::filter_SNPratio(vcf_df_SNPratio, min.SNPratio, min.DP, max.DP)
  
  #Create list of chromosome IDs in the way they appear in the VCF file
  chrList <- BSAvis::extract_chrIDs(vcf.list$meta)
  
  #Plot SNP-ratio across the positions of a given chromosome
  BSAvis::shinyPlot_SNPratio(vcf_df_SNPratio_filt, chrList, chrID, chr, min.SNPratio, degree, span, ranges)
}

