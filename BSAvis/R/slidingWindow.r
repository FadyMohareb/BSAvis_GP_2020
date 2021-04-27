#Calculate sliding windows and mean values (M/WT bulks)

#' @title Sliding Windows
#' @description This function allows to calculate the sliding window based on the chromosome length of the desired chromosome to analyse. 
#' A data frame containing the mean SNP-index values of both bulks gets generated.
#'
#' @param meta meta information stored inside the vcf file
#' @param chrList list of chromosome IDs
#' @param chrID chromosome ID of interest
#' @param windowSize window size (default=1000000)
#' @param windowStep window step (default=10000)
#' @param vcf.df.SNPindex.filt filtered SNP-index data frame
#'
#' @return SNP-index windows
#'
#' @importFrom dplyr %>%
#' @export slidingWindow 
#' @examples
#' ## Default parameters
#' SNPindex_windows <- slidingWindow(meta=vcf_list$meta, 
#'                                   chrList=chromList, 
#'                                   chrID="SL4.0ch03",  
#'                                   vcf.df.SNPindex.filt=vcf_df_SNPindex_filt)
#' ## Custom parameters
#' SNPindex_windows <- slidingWindow(meta=vcf_list$meta, 
#'                                   chrList=chromList, 
#'                                   chrID="SL4.0ch03", 
#'                                   windowSize=2000000, 
#'                                   windowStep=20000, 
#'                                   vcf.df.SNPindex.filt=vcf_df_SNPindex_filt)


slidingWindow <- function(meta, chrList, chrID, windowSize=1000000, windowStep=10000, vcf.df.SNPindex.filt){
  
  # -------------------------- Length of chosen chromosome --------------------------- 
  #Extract those meta lines containing the length of the chromosomes
  lengthLines <- meta[grep("<ID.*length=", meta)]
  #Remove characters before length number (by replacing with "")
  chrLengths <- sub(".*length=", "", lengthLines)
  #Remove characters after length number (by replacing with "")
  chrLengths <- as.numeric(sub(">", "", chrLengths))
  
  #Find length of chosen chromosome (index corresponding to chosen chromosome is the same in chrLengths and chrList)
  chrLength <- chrLengths[grep(chrID, chrList)]
  
  # -------------------- Find start and stop points of each window -------------------- 
  #Find the start points of each window
  windowStart <- seq(from = 1, to = chrLength, by = windowStep)
  #Add window size to each start point to find stop points
  windowStop <- windowStart + windowSize
  
  #Remove windows whose stop positions fall past the chromosome length 
  windowStart <- windowStart[which(windowStop < chrLength)]
  windowStop <- windowStop[which(windowStop < chrLength)]
  
  #Store in data frame start, stop and mid positions for each window
  SNPindex.windows <- data.frame(start = windowStart, stop = windowStop, 
                                 mid = windowStart + (windowStop-windowStart)/2)
  
  # ---------------- Apply sliding window to calculate mean SNP-index -----------------
  #Add new columns to store mean SNP-index (both wild-type and mutant) relative to each window
  SNPindex.windows$mean_SNPindex.WT <- NA
  SNPindex.windows$mean_SNPindex.M <- NA
  
  #Filter data frame by chosen chromosome
  if (length(grep(chrID, chrList)) != 0) { #True if a correct chromosome ID is entered
    df.chr <- vcf.df.SNPindex.filt %>% dplyr::filter(ChromKey==grep(chrID, chrList))
  }
  else {
    #Stop program if an incorrect chromosome ID was entered
    stop("The entered 'chrID' was not found in the VCF file. Please, enter a correct one.")
  }
  
  
  for (n in 1:nrow(SNPindex.windows)) {
    #Restrict data frame to rows whose positions are between the start and stop of the data frame
    df.window <- df.chr[which(df.chr$POS >= SNPindex.windows$start[n] & df.chr$POS <= SNPindex.windows$stop[n]),]
    
    #Calculate mean SNP-index of the variants in that window
    mean_SNPindex.WT <- mean(df.window$SNPindex.WT)
    mean_SNPindex.M <- mean(df.window$SNPindex.M)
    
    #Replace NaN in case they are introduced (when no variants are found in the window)
    if (is.nan(mean_SNPindex.WT)) {mean_SNPindex.WT <- 0.5}
    if (is.nan(mean_SNPindex.M)) {mean_SNPindex.M <- 0.5}
    
    #Set corresponding mean SNP-index of each row
    SNPindex.windows$mean_SNPindex.WT[n] <- mean_SNPindex.WT
    SNPindex.windows$mean_SNPindex.M[n] <- mean_SNPindex.M
  }
  
  return(SNPindex.windows)
}