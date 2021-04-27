#Plot SNP-ratio

#' @title Plot SNP-ratio
#' @description This function allows plotting the SNP-ratio (\emph{Wachsman et al., 2017}).
#' The SNP-ratio gets filtered by the minimum SNP-ratio value (set by the user), before performing LOESS smoothing.
#' The plot gets generated based on the chromosome location of each SNP against the LOESS-fitted SNP-ratio.
#' Whenever the DPI parameter gets specified, the plot will be saved in .TIF format with the specified resolution (dpi value).
#' 
#' \strong{Note:} \emph{for journal publications the preferred format is .TIF, and the minimum resolution is of 600 dpi or 1200 dpi. 
#' Always refer to journal-specific guidelines.}
#' 
#' @references Wachsman et al., 2017
#' 
#' @param vcf.df.SNPratio.filt filtered data frame (containing both bulks)
#' @param chrList list of chromosome IDs
#' @param chrID chromosome ID of interest
#' @param chr chromosome number to print on plot
#' @param min.SNPratio min SNP ratio threshold
#' @param degree LOESS smoothing degree (default=2)
#' @param span LOESS smoothing span (default=0.07)
#' @param filename file name under which the file will be saved (default="plot_SNPratio_chX")
#' @param path path where the file will be saved (default=current working directory)
#' @param dpi resolution value. If no value is given, plots will be generated but not saved
#' @param width width value (default=7.5)
#' @param height height value (default=5)
#' @param units size units (default="in")
#'
#' @importFrom dplyr %>%
#' @export plot_SNPratio 
#' 
#' @examples
#' ## Use default values WITHOUT saving the plot
#' plot_SNPratio(vcf.df.SNPratio.filt=vcf_df_SNPratio_filt, 
#'               chrList=chromList, 
#'               chrID="SL4.0ch03", 
#'               chr=3, 
#'               min.SNPratio=0.3) 
#' ## OR use default values AND save the plot
#' plot_SNPratio(vcf.df.SNPratio.filt=vcf_df_SNPratio_filt, 
#'               chrList=chromList, 
#'               chrID="SL4.0ch03", 
#'               chr=3, 
#'               min.SNPratio=0.3,
#'               dpi=1200) 
#' ## OR customise default parameters
#' plot_SNPratio(vcf.df.SNPratio.filt=vcf_df_SNPratio_filt, 
#'               chrList=chromList, 
#'               chrID="SL4.0ch03", 
#'               chr=3, 
#'               min.SNPratio=0.3, 
#'               degree=2, 
#'               span=0.3, 
#'               filename="SNPratio_ch3", 
#'               path="Document/Plots", 
#'               dpi=1200, 
#'               width=20, 
#'               height=12, 
#'               units="cm")


plot_SNPratio <- function(vcf.df.SNPratio.filt, chrList, chrID, chr, min.SNPratio, degree=2, span=0.07, filename=paste0("plot_SNPratio_ch", chr), path=getwd(), dpi, width=7.5, height=5, units="in") {
  
  #If it is the case, print messages to let user know that default values are being used 
  if (degree==2 & span==0.07) {
    message("=> Applying LOESS smoothing with DEFAULT degree and span values (Degree = 2, Span = 0.07)...")
  }
  else if (degree==2) {
    message("=> Applying LOESS smoothing with DEFAULT degree value (Degree = 2)...")
  }
  else if (span==0.07) {
    message("=> Applying LOESS smoothing with DEFAULT span value (Span = 0.07)...")
  }
  
  #Filter by given chromosome
  if (length(grep(chrID, chrList)) != 0) { #True if a correct chromosome ID is entered
    vcf.df.SNPratio.ch <- vcf.df.SNPratio.filt %>% dplyr::filter(ChromKey==grep(chrID, chrList))
  }
  else {
    #Stop program if an incorrect chromosome ID is entered
    stop("The entered 'chrID' was not found in the VCF file. Please, enter a correct one.")
  }
  
  #Apply LOESS smoothing to SNP-ratio
  loess.fitted <- stats::loess(SNPratio ~ POS, degree=degree, span=span, data=vcf.df.SNPratio.ch)
  
  #Plot the chromosomal location of each SNP against the LOESS-fitted SIMPLE-ratio
  p <- vcf.df.SNPratio.ch %>% dplyr::mutate(smooth = loess.fitted$fitted) %>%
    ggplot2::ggplot(aes(POS, smooth)) +
    ggplot2::geom_point(size = 0.3) +
    ggplot2::geom_hline(yintercept=min.SNPratio, lty=2) +
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("SNP-ratio") +
    ggplot2::labs(title = paste("Chromosome", chr)) 
  
  if(missing(dpi)) {
    #Print message to let user know that the plot was not saved and show required arguments to save it
    message("SNP-ratio plot is being displayed. In order to save it, please specify dpi (and height and width, if different values from the default ones are desired).")
    
    #Show plot
    p
  }
  
  else {
    #Show plot
    p
    
    #Save plot
    ggplot2::ggsave(filename = paste0(filename, ".tiff"), path = path, device = "tiff", dpi = dpi, width = width, height = height, units = units)
    
    #Print messages
    if (width==7.5 & height==5) {
      message("Plot was saved with DEFAULT width and height values (width = 7.5 inches, height = 5 inches).") 
    }
    else if (width==7.5) {
      message("Plot was saved with DEFAULT width value (width = 7.5 inches).")  
    }
    else if (height==5) {
      message("Plot was saved with DEFAULT height value (height = 5 inches).")  
    }
    
    if (path==getwd()) {
      message("Plot was saved in current working directory.")
    }
  }
}