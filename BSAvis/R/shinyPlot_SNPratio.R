#Shiny: SNP-ratio

#' @title Shiny: Plot SNP-ratio 
#' @description Please note that this function is not to be run manually.
#' Creates plots specific to the R-shiny application, showing SNP-ratio values against the position of a specific chromosome.
#' 
#' @param vcf.df.SNPratio.filt filtered data frame (containing both bulks)
#' @param chrList list of chromosome IDs
#' @param chrID chromosome ID of interest
#' @param chr chromosome number to print on plot
#' @param min.SNPratio min SNP ratio threshold
#' @param degree LOESS smoothing degree (default=2)
#' @param span LOESS smoothing span (default=0.07)
#' @param ranges axes ranges (x,y)
#'
#' @details This function is a variant of the plot_SNPratio() function created to meet the needs of the R-shiny application.
#' The difference with the original is that it includes a "ranges" argument that allows setting the limits of the x and y-axis, to show a zoomed area of the SNP-ratio plot when selected. 
#' If no area is selected, the x and y axis limits will be set to null, showing the entire plot.
#'
#' This function does not include any arguments related to plot saving since this functionality is linked to a saving button found inside the BSAvis R-Shiny application.
#'
#' @importFrom dplyr %>%
#' @export shinyPlot_SNPratio 

shinyPlot_SNPratio <- function(vcf.df.SNPratio.filt, chrList, chrID, chr, min.SNPratio, degree=2, span=0.07, ranges) {
  
  #The argument 'ranges' will allow to zoom in the plots while using the shiny App.
  
  #Filter by given chromosome
  vcf.df.SNPratio.ch <- vcf.df.SNPratio.filt %>% dplyr::filter(ChromKey==grep(chrID, chrList))
  
  #Apply LOESS smoothing to SNP-ratio
  loess.fitted <- stats::loess(SNPratio ~ POS, degree=degree, span=span, data=vcf.df.SNPratio.ch)
  
  #Plot the chromosomal location of each SNP against the LOESS-fitted SNP-ratio
  vcf.df.SNPratio.ch %>% dplyr::mutate(smooth = loess.fitted$fitted) %>%
    ggplot2::ggplot(aes(POS, smooth)) +
    ggplot2::geom_point(size = 0.3) +
    ggplot2::geom_hline(yintercept=min.SNPratio, lty=2) +
    ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("SNP-ratio") +
    ggplot2::labs(title = paste("Chromosome", chr)) 
}