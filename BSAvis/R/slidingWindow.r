#Calculate sliding windows and mean values (M/WT bulks)
#'
#' This function allows to calculate the sliding window and the mean values of both bulks.
#'
#' @param Chrom chromosome ID
#' @param chromLength chromosome length
#' @param windowSize window size
#' @param windowStep window step
#' @param df data frame
#'
#' @return windows
#'
#' @export
#' @examples
#' slidingwindow()


slidingWindow <- function(Chrom, chromLength, windowSize, windowStep, df){
  
  #Find the start points of each window
  windowStart <- seq(from = 1, to = chromLength, by = windowStep)
  #Add window size to each start point 
  windowStop <- windowStart + windowSize
  
  #Remove windows whose stop positions fall past the chromosome length 
  windowStart <- windowStart[which(windowStop < chromLength)]
  windowStop <- windowStop[which(windowStop < chromLength)]
  
  #Store in dataframe start, stop and mid positions for each window
  windows <- data.frame(start = windowStart, stop = windowStop, 
                        mid = windowStart + (windowStop-windowStart)/2)
  
  #Add new columns to store mean SNP-index (both wild-type and mutant) relative to each window
  windows$mean_SNPindex.WT <- NA
  windows$mean_SNPindex.M <- NA
  
  #Filter dataframe by chromosome
  df.chrom <- df %>% filter(ChromKey==(Chrom+1))
  
  for (n in 1:nrow(windows)) {
    #Restrict dataframe to rows whose positions are between the start and stop of the dataframe
    df.window <- df.chrom[which(df.chrom$POS >= windows$start[n] & df.chrom$POS <= windows$stop[n]),]
    
    #Calculate mean SNPindex of the variants in that window
    mean_SNPindex.WT <- mean(df.window$SNPindex.WT)
    mean_SNPindex.M <- mean(df.window$SNPindex.M)
    
    #Replace NaN in case they are introduced (when no variants are found in the window)
    if (is.nan(mean_SNPindex.WT)) {mean_SNPindex.WT <- 0.5}
    if (is.nan(mean_SNPindex.M)) {mean_SNPindex.M <- 0.5}
    
    #Set corresponding mean SNP-index of each row
    windows$mean_SNPindex.WT[n] <- mean_SNPindex.WT
    windows$mean_SNPindex.M[n] <- mean_SNPindex.M
  }
  
  return(windows)
}
