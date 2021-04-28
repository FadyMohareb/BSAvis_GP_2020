#Plot delta(SNP Index)
#' @title Plot delta(SNP index)
#' @description This function allows plotting the delta(SNP-index) values against the mid position of the corresponding window of a specific chromosome. 
#' By setting the dpi parameter (resolution), the plot will be automatically saved in .TIF format.
#' 
#' \strong{Note:} \emph{for journal publications the preferred format is .TIF, and the minimum resolution is of 600 dpi or 1200 dpi. 
#' Always refer to journal-specific guidelines.}
#'
#' @param deltaSNPindex.windows filtered data frame (containing both bulks)
#' @param chr chrosome ID
#' @param filename file name under which the file will be saved (default="plot_deltaSNPindex_chX")
#' @param path path where the file will be saved (default=current working directory)
#' @param dpi resolution value. If no value is given, plots will be generated but not saved
#' @param width width value (default=7.5)
#' @param height height value (default=5)
#' @param units size units (default="in")
#'
#' @details ∆(SNP-index) values are plotted against window mid position in a ggplot2line plot. If the dpi argument is not passed, the plot will be shown in the plot panel; however, if a value is added to the dpi argument, the graph will be saved in TIFF format. 
#' Additionally, the name with which to save the file, the directory, the height and width of the plot and their units can be specified in the arguments. 
#' The plot will be saved with default values if different ones are not specified.
#'
#' @export plot_deltaSNPindex
#' @examples
#' ## Use default values WITHOUT saving the plot
#' plot_deltaSNPindex(deltaSNPindex.windows=deltaSNPindex_windows, 
#'                    chr=3)
#' ## OR use default values AND save the plot
#' plot_deltaSNPindex(deltaSNPindex.windows=deltaSNPindex_windows, 
#'                    chr=3,
#'                    dpi=1200)
#' ## OR customise default parameters
#' plot_deltaSNPindex(deltaSNPindex.windows=deltaSNPindex_windows, 
#'                    chr=3,
#'                    filename="deltaSNPindex_ch3", 
#'                    path="Document/Plots", 
#'                    dpi=1200, 
#'                    width=20, 
#'                    height=12, 
#'                    units="cm")


plot_deltaSNPindex <- function(deltaSNPindex.windows, chr, filename=paste0("plot_deltaSNPindex_ch", chr), path=getwd(), dpi, width=7.5, height=5, units="in"){
  
  #Plot delta(SNP-index)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data=deltaSNPindex.windows, aes(mid, delta_SNPindex), size=0.75) +  
    ggplot2::geom_hline(yintercept = 0, lty=2) + 
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(-1, 1) + 
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("delta(SNP-index)") +
    ggplot2::labs(title = paste("Chromosome", chr)) 
  
  if(missing(dpi)) {
    #Print message to let user know that the plot was not saved and show required arguments to save it
    message("Delta(SNP-index) plot is being displayed. In order to save it, please specify dpi (and height and width, if different values from the default ones are desired).")
    
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