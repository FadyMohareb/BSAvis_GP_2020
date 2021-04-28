#SNP-index: wrapper function

#' @title SNP-index Wrapper Function
#' @description This wrapper function is used to fully run the SNP-index method, by calling all the funcitons involved in plotting the mean SNP-index values from each bulk against the mid position of the corresponding window of a specific chromosome. 
#'
#' \deqn{SNPindex=AD_alt/(AD_ref + AD_alt)} 
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
#' @param filename file name under which the file will be saved (default="plot_SNPindex_chX")
#' @param path path where the file will be saved (default=current working directory)
#' @param dpi resolution value. If no value is given, plots will be generated but not saved
#' @param width width value (default=7.5)
#' @param height height value (default=5)
#' @param units size units (default="in")
#' 
#' @details Wrapper function that sequentially calls the required functions involved in generating the SNP-index plot:
#' calc_SNPindex(), filter_SNPindex(), extract_chrIDs(), slidingWindow() and plot_SNPindex(). 
#' The resulting plot will show (for each bulk) mean SNP-index values against the mid position of the corresponding window of a specific chromosome.
#'
#' @importFrom dplyr %>%
#' @export SNPindex_plot
#' @examples
#' ## Use default values WITHOUT saving the plot
#' SNPindex_plot(vcf.list=vcf_list, 
#'          wtBulk="pool_S3781_minus", 
#'          mBulk="pool_S3781_plus", 
#'          chrID="SL4.0ch03", 
#'          chr=3)
#' ## OR use default values AND save the plot
#' SNPindex_plot(vcf.list=vcf_list, 
#'          wtBulk="pool_S3781_minus", 
#'          mBulk="pool_S3781_plus", 
#'          chrID="SL4.0ch03", 
#'          chr=3,
#'          dpi=1200)
#' #OR customise default parameters
#' SNPindex_plot(vcf.list=vcf_list, 
#'          wtBulk="pool_S3781_minus", 
#'          mBulk="pool_S3781_plus", 
#'          variants="all",
#'          min.SNPindex=0.25, 
#'          max.SNPindex=0.8, 
#'          min.DP=60, 
#'          max.DP=250, 
#'          min.GQ=98,
#'          chrID="SL4.0ch03", 
#'          chr=3,
#'          windowSize=2000000, 
#'          windowStep=20000,
#'          filename="SNPindex_chrom03", 
#'          path="Document/Plots", 
#'          dpi=1200,
#'          width=20, 
#'          height=12, 
#'          units="cm")


SNPindex_plot <- function(vcf.list, wtBulk, mBulk, variants="SNP",
                     min.SNPindex=0.3, max.SNPindex=0.9, min.DP=50, max.DP=200, min.GQ=99,
                     chrID, chr, windowSize=1000000, windowStep=10000, 
                     filename=paste0("plot_SNPindex_ch", chr), path=getwd(), 
                     dpi, width=7.5, height=5, units="in"){
  
  #Calculate SNP-index of each variant in each bulk
  vcf_df_SNPindex <- BSAvis::calc_SNPindex(vcf.list$df, wtBulk, mBulk, variants)
  
  #Filter variants
  vcf_df_SNPindex_filt <- BSAvis::filter_SNPindex(vcf_df_SNPindex, min.SNPindex, max.SNPindex, min.DP, max.DP, min.GQ)
  
  #Create list of chromosome IDs in the way they appear in the VCF file
  chrList <- BSAvis::extract_chrIDs(vcf.list$meta)
  
  #Apply sliding window to calculate mean SNP-index in each window of a specifc size and step for each bulk
  SNPindex_windows <- BSAvis::slidingWindow(vcf.list$meta, chrList, chrID, windowSize, windowStep, vcf_df_SNPindex_filt)
  
  #Plot SNP-index across the positions of a given chromosome
  BSAvis::plot_SNPindex(SNPindex_windows, chr, filename, path, dpi, width, height, units)
}

