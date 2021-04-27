#Shiny: Plot delta(SNP Index)
#' @title Plot delta(SNP index)
#' @description Please note that this function is not to be run manually.
#'
#' @param vcf.df.window.delta filtered data frame (containing both bulks)
#' @param chr chrosome ID
#' @param ranges axes ranges (x,y)
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
