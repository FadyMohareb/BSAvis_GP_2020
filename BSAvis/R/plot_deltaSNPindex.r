#Plot delta(SNP Index)
#'
#' This function allows plotting delta(SNP Index).
#'
#' @param vcf.df.bulks.filt filtered data frame (containing both bulks)
#' @param Chrom chromosome ID
#'
#' @export
#' @examples
#' plot_deltaSNPindex()

plot_deltaSNPindex <- function(vcf.df.bulks.filt, Chrom){
  
  #Plot delta SNP index
  ggplot2::ggplot() +
    ggplot2::geom_line(data=vcf.df.bulks.filt, aes(mid, delta_SNPindex), size=0.75) +  
    ggplot2::geom_hline(yintercept = 0, lty=2) + 
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(-1, 1) + 
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("delta(SNP-index)") +
    ggplot2::labs(title = paste("Chromosome", Chrom)) 
}
