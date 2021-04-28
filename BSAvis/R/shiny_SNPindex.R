#Shiny: SNP-index wrapper function

#' @title SNP-index Wrapper Function
#' @description Please note that this function is not to be run manually.
#' Wrapper function for SNP-index plots specific to the R-shiny application. 
#' Calls all functions involved in plotting the mean SNP-index values from both bulks (or just one of them) against the mid position of the corresponding window of a specific chromosome.
#'
#' @param vcf.list object containing meta information and vcf data frame
#' @param wtBulk Wild-Type pool
#' @param mBulk Mutant pool 
#' @param variants variants to be considered. Default is "SNP" (allowed: "SNP" or "all")
#' @param min.SNPindex min value allowed for the SNP index (default=0.3)
#' @param max.SNPindex max value allowed for the SNP index (default=0.9)
#' @param min.DP min value allowed for the read depth (default=50)
#' @param max.DP max value allowed for the read depth (default=200)
#' @param min.GQ min value allowed for the genotype quality (default=99)
#' @param chrID chromosome ID of interest
#' @param chr chromosome name printed on the plot
#' @param windowSize window size (default=1000000)
#' @param windowStep window step (default=10000)
#' @param bulk bulk/bulks to take into consideration
#' @param ranges axes ranges (x,y)
#' 
#' @details Wrapper function for generating the SNP-index plots, inside the R-Shiny application.
#' The involved functions are the following (in oder): calc_SNPindex(), filter_SNPindex(), extract_chrIDs(), slidingWindow() and shinyPlot_SNPindex().
#' 
#' This function is not intended to be run by the user, as it was specifically created to be used by the BSAvis_shiny() function.
#' @importFrom dplyr %>%
#' @export shiny_SNPindex


shiny_SNPindex <- function(vcf.list, wtBulk, mBulk, variants,
                           min.SNPindex, max.SNPindex, min.DP, max.DP, min.GQ,
                           chrID, chr, windowSize, windowStep, bulk, ranges){
 
  #Calculate SNP-index of each variant in each bulk
  vcf_df_SNPindex <- BSAvis::calc_SNPindex(vcf.list$df, wtBulk, mBulk, variants)
 
  #Filter variants
  vcf_df_SNPindex_filt <- BSAvis::filter_SNPindex(vcf_df_SNPindex, min.SNPindex, max.SNPindex, min.DP, max.DP, min.GQ)
 
  #Create list of chromosome IDs in the way they appear in the VCF file
  chrList <- BSAvis::extract_chrIDs(vcf.list$meta)
 
  #Apply sliding window to calculate mean SNP-index in each window of a specifc size and step for each bulk
  SNPindex_windows <- BSAvis::slidingWindow(vcf.list$meta, chrList, chrID, windowSize, windowStep, vcf_df_SNPindex_filt)
 
  #Plot SNP-index across the positions of a given chromosome
  BSAvis::shinyPlot_SNPindex(SNPindex_windows, chr, bulk, ranges)
}