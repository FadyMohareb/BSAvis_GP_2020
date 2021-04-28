#Shiny: Plot delta(SNP Index)
#' @title Plot delta(SNP index)
#' @description Please note that this function is not to be run manually.
#' Creates plots specific to the R-shiny application, showing the delta(SNP-index) values against the mid position of the corresponding window of a specific chromosome.

#' @param vcf.df.window.delta filtered data frame (containing both bulks)
#' @param chr chrosome ID
#' @param ranges axes ranges (x,y)
#'
#' @details This function is a variant of the plot_deltaSNPindex() function created to meet the needs of the R-shiny application.
#' The difference with the original is that it includes a "ranges" argument that allows setting the limits of the x and y-axis, to show a zoomed area of the delta(SNP-index) plot when selected. 
#' If no area is selected, the x and y axis limits will be set to null, showing the entire plot.
#'
#' This function does not include any arguments related to plot saving since this functionality is linked to a saving button found inside the BSAvis R-Shiny application.
#'
#' @export shinyPlot_deltaSNPindex

shinyPlot_deltaSNPindex <- function(vcf.df.window.delta, chr, ranges){
  
  #The argument 'ranges' will allow to zoom in the plots while using the shiny App.
  
  #Plot delta(SNP-index)
  ggplot2::ggplot() +
    ggplot2::geom_line(data=vcf.df.window.delta, aes(mid, delta_SNPindex), size=0.75) +  
    ggplot2::geom_hline(yintercept = 0, lty=2) + 
    ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(-1, 1) + 
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("delta(SNP-index)") +
    ggplot2::labs(title = paste("Chromosome", chr)) 
}
