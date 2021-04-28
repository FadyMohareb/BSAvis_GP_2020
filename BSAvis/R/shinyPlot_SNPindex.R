#Shiny: SNP-indices

#' @title Shiny: Plot SNP-indices
#' @description Please note that this function is not to be run manually.
#' Creates plots specific to the R-shiny application, showing (for each bulk, or just one of them) the mean SNP-index values against the mid position of the corresponding window of a specific chromosome.
#' 
#' @param SNPindex.windows filtered data frame (containing both bulks)
#' @param chr chrosome ID
#' @param bulk bulk/bulks to take into consideration when plotting
#' @param ranges axes ranges (x,y)
#' 
#' @details This function is a variant of the plot_SNPindex() function created to meet the needs of the R-shiny application.
#' Firstly, a ggplot object is created with all the necessary plotting elements, exept for the lines specific for each bulk (corresponding to mean SNP-index against the mid position of the corresponding window of a specific chromosome).
#' Based on the value of the bulk argument, a SNP-index line/s will be added to the previously created ggplotobject. 
#'
#' If bulk=0, two lines corresponding to the SNP-indices of both wild-type and mutant bulks will be added.
#' If bulk=1 or bulk=2, a single line will be added to the plot, respectively corresponding to the wild-type or mutant bulk SNP-indices.
#' 
#' The value of the bulk parameter will be set inside the R-shiny application function (BSAvis_shiny()) depending on the user's interactive-input.
#' Additionally, this function includes a "ranges" argument that allows setting the limits of the x and y-axis, to show a zoomed area of the SNP-index plot when selected. 
#' If no area is selected, the x and y axis limits will be set to null, showing the entire plot.
#'
#' This function does not include any arguments related to plot saving since this functionality is linked to a saving button found inside the BSAvis R-Shiny application.
#'
#' @export shinyPlot_SNPindex

shinyPlot_SNPindex <- function(vcf.df.window.SNPindex, chr, bulk, ranges) {
 
  #This function includes the argument 'bulk' in order to add geom_line() to ggplot()
  #to show just the line/s for the chosen bulk/s by the user in the shiny app.
  #The possible values of this argument and the corresponding meanings are:
  #bulk=0 - visualise SNP-index line from WT and M pools
  #bulk=1 - visualise SNP-index line from WT pool
  #bulk=2 - visualise SNP-index line from M pools
 
  #The argument 'ranges' will allow to zoom in the plots while using the shiny App.
 
  #Plot SNP-indexes of both pools.
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0.5, lty=2) +
    ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("Position (Mb)") +
    ggplot2::ylab("SNP-index") +
    ggplot2::labs(title = paste("Chromosome", chr)) +
    ggplot2::theme(legend.title = element_blank())
 
  if (bulk==0) {
    p <- p + ggplot2::geom_line(data=vcf.df.window.SNPindex, aes(mid, mean_SNPindex.WT, color="WT pool"), size=0.75) +
      ggplot2::geom_line(data=vcf.df.window.SNPindex, aes(mid, mean_SNPindex.M, color="Mutant pool"), size=0.75)
  }
 
  else if (bulk==1) {
    p <- p + ggplot2::geom_line(data=vcf.df.window.SNPindex, aes(mid, mean_SNPindex.WT, color="WT pool"), size=0.75) +
             scale_color_manual(values=c("#00BFC4"))
  }
 
  else if (bulk==2) {
    p <- p + ggplot2::geom_line(data=vcf.df.window.SNPindex, aes(mid, mean_SNPindex.M, color="Mutant pool"), size=0.75) +
             scale_color_manual(values=c("#F8766D"))
  }
 
  #Show plot
  p
}