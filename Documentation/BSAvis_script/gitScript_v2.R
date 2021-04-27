#####################
###R-Script BSAvis###
#####################

#Read vcf file and convert vcfR object to dataframe
readBSA_vcf <- function(file) {
  
  vcf <- vcfR::read.vcfR(file, verbose = FALSE)

  vcf.df <- vcfR::vcfR2tidy(vcf,
                            format_fields = c("AD", "DP", "GQ", "GT"),  
                            gt_column_prepend = "")
  
  #Limit dataframe to columns of interest
  vcf.df <- vcf.df$gt
  
  #Split AD column into reference and alterate AD values
  vcf.df <- vcf.df %>% tidyr::separate(AD, sep = ",", into = c("AD_ref", "AD_alt"))
  
  vcf.list <- list(meta=vcf@meta, df=vcf.df)
  
  return(vcf.list)
}


calc_SNPindex <- function(vcf.df, WT.bulk, M.bulk) {
  
  #Create dataframe for each bulk sample AND include a new column with SNP index
  vcf.df.WTbulk <- vcf.df %>% dplyr::filter(Indiv==WT.bulk) %>% dplyr::mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))
  vcf.df.Mbulk <- vcf.df %>% dplyr::filter(Indiv==M.bulk) %>% dplyr::mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))

  vcf.df.jointBulks <- dplyr::inner_join(vcf.df.WTbulk, 
                                         vcf.df.Mbulk, 
                                         by = c("ChromKey","POS"), 
                                         copy = F, 
                                         suffix = c(".WT", ".M"))
  
  return(vcf.df.jointBulks)
}


filter_variants <- function(vcf.df.bulks, min_SNPindex=0.3, max_SNPindex=0.9, min_DP=50, max_DP=200, min_GQ=99){
  
  #Filter by min SNP-index
  vcf.df.bulks <- vcf.df.bulks %>% dplyr::filter(SNPindex.WT >= min_SNPindex & SNPindex.M >= min_SNPindex &
                                                 SNPindex.WT <= max_SNPindex & SNPindex.M <= max_SNPindex) 
  #Filter by min and max depth
  vcf.df.bulks <- vcf.df.bulks %>% dplyr::filter(DP.WT >= min_DP & DP.M >= min_DP &
                                                 DP.WT <= max_DP & DP.M <= max_DP) 
  #Filter by min genotype quality
  vcf.df.bulks.filt <- vcf.df.bulks %>% dplyr::filter(GQ.WT >= min_GQ & GQ.M >= min_GQ) 
  
  return(vcf.df.bulks.filt)
}


plot_SNPindex <- function(vcf.df.bulks.filt, Chrom) {
  
  #Plot SNP indexes of both pools
  ggplot2::ggplot() +
    ggplot2::geom_line(data=vcf.df.bulks.filt, aes(mid, mean_SNPindex.WT, color="WT pool"), size=0.75) +  
    ggplot2::geom_line(data=vcf.df.bulks.filt, aes(mid, mean_SNPindex.M, color="Mutant pool"), size=0.75) +
    ggplot2::geom_hline(yintercept = 0.5, lty=2) + 
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::scale_x_continuous(labels=function(x)x/1000000) +
    ggplot2::ylim(0, 1) + 
    ggplot2::xlab("Position (Mb)") + 
    ggplot2::ylab("SNP-index") +
    ggplot2::labs(title = paste("Chromosome", Chrom)) + 
    ggplot2::theme(legend.title = element_blank()) 
}


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


extract_chromLength <- function(chromID, meta){
  
  #Extract meta line containing length info on the specific chromosome
  lengthLine <- meta[grep(chromID, meta)]
  #Remove (by replacing with "") previous characters to length number
  length <- sub(".*length=", "", lengthLine)
  #Set chromosome size (and remove ">" character which appears after length number)
  chromLength <- as.numeric(sub(">", "", length))
  
  return(chromLength)
}


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


calc_deltaSNPindex <- function(vcf.df.window) {
  
  #Calculate delta SNP-index and add it as new column
  vcf.df.deltaSNPindex <- vcf.df.window %>% dplyr::mutate("delta_SNPindex" = (mean_SNPindex.M - mean_SNPindex.WT))
  return(vcf.df.deltaSNPindex)
}
