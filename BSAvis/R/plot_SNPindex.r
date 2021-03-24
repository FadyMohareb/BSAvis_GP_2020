#Plot (SNP Index)
#'
#' This function allows plotting the (SNP Index).
#'
#' @param vcf.df.bulks.filt filtered data frame (containing both bulks)
#' @param Chrom chromosome ID
#'
#' @export
#' @examples
#' plot_SNPindex()

plot_SNPindex <- function(vcf.df.bulks.filt, Chrom) {
  
  #Plot SNP indexes of both pools
  ggplot2::ggplot() +
    ggplot2::geom_line(data=vcf.df.bulks, aes(mid, mean_SNPindex.WT, color="WT pool"), size=0.75) +  
    ggplot2::geom_line(data=vcf.df.bulks, aes(mid, mean_SNPindex.M, color="Mutant pool"), size=0.75) +
    ggplot2::geom_hline(yintercept = 0.5, lty=2) + 
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(0, 1) + 
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("SNP-index") +
    ggplot2::labs(title = paste("Chromosome", Chrom)) + 
    ggplot2::theme(legend.title = element_blank()) 
}

