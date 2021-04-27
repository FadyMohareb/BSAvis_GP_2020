#Plot SNP-indices

#' @title Plot SNP-indices
#' @description This function allows plotting the SNP-indices.
#' Whenever the dpi parameter gets specified, the plot will be saved in .TIF format with the specified resolution (dpi value).
#' 
#' \strong{Note:} \emph{for journal publications the preferred format is .TIF, and the minimum resolution is of 600 dpi or 1200 dpi. 
#' Always refer to journal-specific guidelines.}
#' 
#' @param SNPindex.windows filtered data frame (containing both bulks)
#' @param chr chrosome ID
#' @param filename file name under which the file will be saved (default="plot_SNPindex_chX")
#' @param path path where the file will be saved (default=current working directory)
#' @param dpi resolution value. If no value is given, plots will be generated but not saved
#' @param width width value (default=7.5)
#' @param height height value (default=5)
#' @param units size units (default="in")
#' 
#' @export plot_SNPindex 
#' @examples
#' ## OR use default values AND save the plot
#' plot_SNPindex(SNPindex.windows=SNPindex_windows, 
#'               chr=3,
#'               dpi=1200) 
#' ## OR customise default parameters
#' plot_SNPindex(SNPindex.windows=SNPindex_windows, 
#'               chr=3, 
#'               filename="SNPindex_ch3", 
#'               path="Document/Plots", 
#'               dpi=1200, 
#'               width=20, 
#'               height=12, 
#'               units="cm") 


plot_SNPindex <- function(SNPindex.windows, chr, filename=paste0("plot_SNPindex_ch", chr), path=getwd(), dpi, width=7.5, height=5, units="in") {
  
  #Plot SNP-indexes of both pools
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data=SNPindex.windows, aes(mid, mean_SNPindex.WT, color="WT pool"), size=0.75) +  
    ggplot2::geom_line(data=SNPindex.windows, aes(mid, mean_SNPindex.M, color="Mutant pool"), size=0.75) +
    ggplot2::geom_hline(yintercept = 0.5, lty=2) + 
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(0, 1) + 
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("SNP-index") +
    ggplot2::labs(title = paste("chrosome", chr)) + 
    ggplot2::theme(legend.title = element_blank())
  
  if(missing(dpi)){
    #Print message to let user know that the plot was not saved and show required arguments to save it
    message("SNP-index plot is being displayed. In order to save it, please specify dpi (and height and width if different values from the default ones are desired).")
    
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