#Read vcf file

#' @title Read vcf file
#' @description This function allows to read a vcf file, converting it into a data frame. 
#' Format fields get separated as follows: AD, DP, GQ, GT and GT_alleles. 
#' 
#' Allelic depths (AD) get split into reference and alterate AD values ("AD_ref" and "AD_alt").
#' 
#' @param file VCF file containing wild-type and mutant pools as samples, which must include the AD column.
#'
#' @return a list containing a character vector (with meta information) and a data frame (with variant information).
#' @details Using the vcf R package, a VCF file is read and convertedinto a data frame containing the followingc olumns:
#' chromosome, position, individual (sample ID in the VCF file), reference AD, alternate AD, DP, GQ, GT, and GT_alleles. 
#' The meta information of the VCF file is also extracted. 
#' The output is a list containing two elements: meta information and the mentioned dataframe
#''
#' @importFrom dplyr %>%
#' @export readBSA_vcf
#' 
#' @examples
#' ## Read vcf file
#' vcf_list <- readBSA_vcf("dataset1_pools.vcf")


readBSA_vcf <- function(file) {
  
  vcf <- vcfR::read.vcfR(file, verbose = FALSE)

  vcf.df <- vcfR::vcfR2tidy(vcf,
                            format_fields = c("AD", "DP", "GQ", "GT"),  
                            gt_column_prepend = "")
  
  #Limit dataframe to columns of interest
  vcf.df <- vcf.df$gt
  
  #Split AD column into reference and alterate AD values
  vcf.df <- vcf.df %>% tidyr::separate(AD, sep = ",", into = c("AD_ref", "AD_alt"))
  
  #Create list containing meta information and data frame
  vcf.list <- list(meta=vcf@meta, df=vcf.df)
  
  return(vcf.list)
}
