#SNP-ratio: wrapper function

#' @title SNP-ratio Wrapper Function
#' @description This wrapper function is used to fully run the SNP-ratio method, by calling all the funcitons involved in plotting the SNP-ratio values across positions of a specific chromosome.
#'
#' \deqn{SNPindex=AD_alt/(AD_ref + AD_alt)} 
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
#' @param filename file name under which the file will be saved (default="plot_SNPratio_chX")
#' @param path path where the file will be saved (default=current working directory)
#' @param dpi resolution value. If no value is given, plots will be generated but not saved
#' @param width width value (default=7.5)
#' @param height height value (default=5)
#' @param units size units (default="in")
#'
#' @details Wrapper function that sequentially calls the required functions involved in generating the SNP-ratio plot:
#' calc_SNPratio(), filter_SNPratio(), extract_chrIDs() and plot_SNPratio().
#' The resulting plot will show SNP-ratio values (for both bulks) across chromosome positions.
#'
#' 
#' @importFrom dplyr %>%
#' @export SNPratio_plot
#' @examples
#' ## Use default values WITHOUT saving the plot
#' SNPratio_plot(vcf.list=vcf_list, 
#'          wtBulk="pool_S3781_minus", 
#'          mBulk="pool_S3781_plus", 
#'          chrID="SL4.0ch03", 
#'          chr=3)
#' ## OR use default values AND save the plot
#' SNPratio_plot(vcf.list=vcf_list, 
#'          wtBulk="pool_S3781_minus", 
#'          mBulk="pool_S3781_plus", 
#'          chrID="SL4.0ch03", 
#'          chr=3, 
#'          dpi=1200)
#' ## OR customise default parameters
#' SNPratio_plot(vcf.list=vcf_list, 
#'          wtBulk="pool_S3781_minus", 
#'          mBulk="pool_S3781_plus", 
#'          variants="all",
#'          min.SNPratio=0.3, 
#'          min.DP=60, 
#'          max.DP=250,
#'          chrID="SL4.0ch03", 
#'          chr=3, 
#'          degree=1, 
#'          span=0.05, 
#'          filename="deltaSNPindex_chrom03", 
#'          path="Document/Plots", 
#'          dpi=1200, 
#'          width=20, 
#'          height=12, 
#'          units="cm")


SNPratio_plot <- function(vcf.list, wtBulk, mBulk, variants="SNP",
                     min.SNPratio=0.1, min.DP=50, max.DP=200,
                     chrID, chr, degree=2, span=0.07, 
                     filename=paste0("plot_SNPratio_ch", chr), path=getwd(), 
                     dpi, width=7.5, height=5, units="in"){
  
  #Calculate SNP-ratio of each variant
  vcf_df_SNPratio <- BSAvis::calc_SNPratio(vcf.list$df, wtBulk, mBulk, variants)
  
  #Filter variants
  vcf_df_SNPratio_filt <- BSAvis::filter_SNPratio(vcf_df_SNPratio, min.SNPratio, min.DP, max.DP)
  
  #Create list of chromosome IDs in the way they appear in the VCF file
  chrList <- BSAvis::extract_chrIDs(vcf.list$meta)
  
  #Plot SNP-ratio across the positions of a given chromosome
  BSAvis::plot_SNPratio(vcf_df_SNPratio_filt, chrList, chrID, chr, min.SNPratio, degree, span, filename, path, dpi, width, height, units)
}
