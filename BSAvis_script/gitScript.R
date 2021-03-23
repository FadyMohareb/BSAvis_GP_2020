#####################
###R-Script BSAvis###
#####################

#Load additional packages
library(vcfR)
library(tidyr)
library(dplyr)
library(ggplot2)

#Read vcf file and convert vcfR object to dataframe
readBSA_vcf <- function(file) {
  
  vcf <- vcfR::read.vcfR(file, verbose = FALSE)

  vcf_tidy <- vcfR2tidy(vcf,
                        format_fields = c("AD", "DP", "GQ", "GT"),  
                        gt_column_prepend = "")
  
  #Limit dataframe to columns of interest
  vcf_tidy <- vcf_tidy$gt
  
  #Split AD column into reference and alterate AD values
  vcf_tidy <- vcf_tidy %>% separate(AD, sep = ",", into = c("AD_ref", "AD_alt"))
  return(vcf_tidy)
}


calc_SNPindex <- function(file, WTbulk, Mbulk) {
  
  #Call function to read vcf file and obtain dataframe
  vcf.df <- readBSA_vcf(file)
  
  #Create dataframe for each bulk sample AND include a new column with SNP index
  vcf.df.WTbulk <- vcf_tidy %>% filter(Indiv==WTbulk) %>% mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))
  vcf.df.Mbulk <- vcf_tidy %>% filter(Indiv==Mbulk) %>% mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))

  vcf.df.jointBulks <- inner_join(vcf.df.WTbulk, 
                                  vcf.df.Mbulk, 
                                  by = c("ChromKey","POS"), 
                                  copy = F, 
                                  suffix = c(".WT", ".M"))
  return(vcf.df.jointBulks)
}


filter_variants <- function(vcf.df.bulks, min_SNPindex=0.3, max_SNPindex=0.9, min_DP=50, max_DP=200, min_GQ=99){
  
  #Filter by min SNP-index
  vcf.df.bulks <- vcf.df.bulks %>% filter(SNPindex.WT >= min_SNPindex & SNPindex.M >= min_SNPindex &
                                          SNPindex.WT <= max_SNPindex & SNPindex.M <= max_SNPindex) 
  #Filter by min and max depth
  vcf.df.bulks <- vcf.df.bulks %>% filter(DP.WT >= min_DP & DP.M >= min_DP &
                                          DP.WT <= max_DP & DP.M <= max_DP) 
  #Filter by min genotype quality
  vcf.df.bulks.filt <- vcf.df.bulks %>% filter(GQ.WT >= min_GQ & GQ.M >= min_GQ) 
  
  return(vcf.df.bulks.filt)
}


plot_SNPindex <- function(vcf.df.bulks, Chrom) {
  
  #Plot SNP indexes of both pools
  ggplot() +
    geom_line(data=vcf.df.bulks, aes(mid, mean_SNPindex.WT, color="WT pool"), size=0.75) +  
    geom_line(data=vcf.df.bulks, aes(mid, mean_SNPindex.M, color="Mutant pool"), size=0.75) +
    geom_hline(yintercept = 0.5, lty=2) + 
    coord_cartesian(expand=FALSE) +
    scale_x_continuous(labels=function(x)x/1000000) +
    ylim(0, 1) + 
    xlab("Position (Mb)") + 
    ylab("SNP-index") +
    labs(title = paste("Chromosome", Chrom)) + 
    theme(legend.title = element_blank()) 
}


plot_deltaSNPindex <- function(vcf.df.bulks, Chrom){
  
  #Plot delta SNP index
  ggplot() +
    geom_line(data=vcf.df.bulks, aes(mid, delta_SNPindex), size=0.75) +  
    geom_hline(yintercept = 0, lty=2) + 
    coord_cartesian(expand=FALSE) +
    scale_x_continuous(labels=function(x)x/1000000) +
    ylim(-1, 1) + 
    xlab("Position (Mb)") + 
    ylab("delta(SNP-index)") +
    labs(title = paste("Chromosome", Chrom)) 
}


extract_chromLength <- function(chromID, vcf){
  
  #Extract meta line containing length info on the specific chromosome
  lengthLine <- vcfMeta[grep(chromID, vcf@meta)]
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
