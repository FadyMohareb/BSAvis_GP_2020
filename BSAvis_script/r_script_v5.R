#####################
###R-Script BSAvis###
#####################

#Clear workspace
rm(list=ls())

#Close any open graphics devices
graphics.off()

#Load additional packages
#install.packages("vcfR")
library(vcfR)
library(tidyr)
#library(tidyverse)
library(dplyr)
library(ggplot2)

#===================================================================================#

#Read and store input vcf file
vcf <- read.vcfR("dataset1.vcf", verbose = FALSE)

#Convert vcfR object to tidy dataframe, suiatble for analysis with dyplr, tidyr, ggplot2
vcf_tidy <- vcfR2tidy(vcf,
                      format_fields = c("AD", "DP", "GQ", "GT"),  
                      gt_column_prepend = "")

#Save vcf_tidy in a new variable so that we have two copies of it. 
#In case we transformed one in a wrong one, we have other copy and 
#avoid re-running the last command which takes long
vcf_tidy_copy <- vcf_tidy
#vcf_tidy <- vcf_tidy_copy

#Store meta info
meta <- vcf@meta

#Limit dataframe to columns of interest
vcf_tidy <- vcf_tidy$gt

#Split AD column into reference and alterate AD values
vcf_tidy <- vcf_tidy %>% separate(AD, sep = ",", into = c("AD_ref", "AD_alt"))

#Create dataframe for each bulk sample AND includE a new column with SNP index
vcf.pool_S3781_minus <- vcf_tidy %>% filter(Indiv=="pool_S3781_minus") %>% mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))
vcf.pool_S3781_plus <- vcf_tidy %>% filter(Indiv=="pool_S3781_plus") %>% mutate("SNPindex"= as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)))


################################## start of FUNCTIONS DEFINITIONS ##################################
#(FUNCTION 1) Define filter function:
filter_variants <- function(poolVCF_df, poolName, min_SNPindex, max_SNPindex, min_DP, max_DP, min_GQ){
  
  poolVCF_df <- poolVCF_df %>% filter(SNPindex > min_SNPindex & 
                                        SNPindex <= max_SNPindex) #filter by min SNP-index
  poolVCF_df <- poolVCF_df %>% filter(DP > min_DP & 
                                        DP < max_DP ) #filter by min and max depth
  poolVCF_df <- poolVCF_df %>% filter(GQ >= min_GQ) #filter by min genotype quality
}


#(FUNCTION 2) Define plot function for SNP-index:
plot_SNPindex <- function(poolsVCF_df, Chrom) {

  #Plot SNP indexes of both pools
  ggplot() +
  geom_line(data=poolsVCF_df, aes(mid, mean_SNPindex.WT, color="WT pool S3781"), size=0.75) +  
  geom_line(data=poolsVCF_df, aes(mid, mean_SNPindex.M, color="Mutant pool S3781"), size=0.75) +
  geom_hline(yintercept = 0.5, lty=2) + 
  coord_cartesian(expand=FALSE) +
  scale_x_continuous(labels=function(x)x/1000000) +
  ylim(0, 1) + 
  xlab("Position (Mb)") + 
  ylab("SNP-index") +
  labs(title = paste("Chromosome", Chrom)) + 
  theme(legend.title = element_blank()) 
  #ggsave(paste0("SNPindex_ch", Chrom, ".png"), width=89, height=50, units="mm")
}


#(FUNCTION 3) Define plot function for delta SNP-index:
plot_deltaSNPindex <- function(poolsVCF_df, Chrom){

  #Plot delta SNP index
  ggplot() +
    geom_line(data=poolsVCF_df, aes(mid, delta_SNPindex), size=0.75) +  
    geom_hline(yintercept = 0, lty=2) + 
    coord_cartesian(expand=FALSE) +
    scale_x_continuous(labels=function(x)x/1000000) +
    ylim(-1, 1) + 
    xlab("Position (Mb)") + 
    ylab("delta(SNP-index)") +
    labs(title = paste("Chromosome", Chrom)) 
    #ggsave(paste0("plot_deltaSNPindex_ch"), Chrom, "_ds1.png"), width=89, height=50, units="mm")
}

#(FUNCTION 4) Define function for extracting chromosome length:
extract_chromLength <- function(chromID, vcfMeta){
  
  #Extract meta line containing length info on the specific chromosome
  lengthLine <- vcfMeta[grep(chromID, vcfMeta)]
  #Remove (by replacing with "") previous characters to length number
  length <- sub(".*length=", "", lengthLine)
  # set chromosome size (and remove ">" character which appears after length number)
  chromLength <- as.numeric(sub(">", "", length))
  return(chromLength)
  
}

#(FUNCTION 5) Define sliding window function:
slidingWindow <- function(chrom, chromLength, windowSize, windowStep, df){
  
  #Find the start points of each window
  windowStart <- seq(from = 1, to = chromLength, by = windowStep)
  #Add window size to each start point 
  windowStop <- windowStart + windowSize
  
  # no windows start before the end of chromosome 8
  sum(windowStart > chromLength)
  # but some window stop positions do occur past the final point
  sum(windowStop > chromLength)
  
  #Remove windows whose stop positions fall past the chromosome length 
  windowStart <- windowStart[which(windowStop < chromLength)]
  windowStop <- windowStop[which(windowStop < chromLength)]
  
  #Store in dataframe start, stop and mid positions for each window
  windows <- data.frame(start = windowStart, stop = windowStop, 
                        mid = windowStart + (windowStop-windowStart)/2)
  
  #Add new columns to store mean SNP-index (both wild-type and mutant) relative to each window
  windows$mean_SNPindex.WT <- NA
  windows$mean_SNPindex.M <- NA
  
  df.chrom <- df %>% filter(ChromKey==(chrom+1))
  
  for (n in 1:nrow(windows)) {
    df.window <- df.chrom[which(df.chrom$POS >= windows$start[n] & df.chrom$POS <= windows$stop[n]),]
    mean_wt_SNPindex <- mean(df.window$SNPindex.S3781_minus)
    mean_m_SNPindex <- mean(df.window$SNPindex.S3781_plus)
    windows$mean_SNPindex.WT[n] <- mean_wt_SNPindex
    windows$mean_SNPindex.M[n] <- mean_m_SNPindex
  }
  return(windows)
  
}
################################### end of FUNCTIONS DEFINITIONS ###################################


#Apply filter function
vcf.pool_S3781_minus.filt <- filter_variants(poolVCF_df=vcf.pool_S3781_minus, 
                                             min_SNPindex=0.3, max_SNPindex=0.9, 
                                             min_DP=50, max_DP=200, min_GQ=99)

vcf.pool_S3781_plus.filt <- filter_variants(poolVCF_df=vcf.pool_S3781_plus, 
                                            min_SNPindex=0.3, max_SNPindex=0.9, 
                                            min_DP=50, max_DP=200, min_GQ=99)

#Join wild-type and mutant dataframes
vcf.joint_pools_S3781 <- inner_join(vcf.pool_S3781_minus.filt, 
                                    vcf.pool_S3781_plus.filt, 
                                    by = c("ChromKey","POS"), 
                                    copy = F, 
                                    suffix = c(".S3781_minus", ".S3781_plus"))

#Extract chromosome length of a specific chromosome
chromLength <- extract_chromLength("SL4.0ch01", vcf@meta)

#Create dataframe with SNP-indexes in each window
pools_S3781_window.df <- slidingWindow(1, chromLength, 2000000, 10000, vcf.joint_pools_S3781)

#Plot SNP-index 
plot_SNPindex(poolsVCF_df=pools_S3781_window.df, Chrom=1) 

#Calculate delta SNP-index and add it as new column
pools_S3781_window.df <- pools_S3781_window.df %>% mutate("delta_SNPindex" = (mean_SNPindex.M - mean_SNPindex.WT))

#Plot delta SNP-index
plot_deltaSNPindex(poolsVCF_df=pools_S3781_window.df, Chrom=1)


#Save image
save.image("r_script_v5.RData")