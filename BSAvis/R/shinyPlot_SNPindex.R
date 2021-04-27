#Shiny: SNP-indices

#' @title Shiny: Plot SNP-indices
#' @description Please note that this function is not to be run manually.
#' 
#' @param SNPindex.windows filtered data frame (containing both bulks)
#' @param chr chrosome ID
#' @param bulk bulk/bulks to take into consideration when plotting
#' @param ranges axes ranges (x,y)
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