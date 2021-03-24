#Read file and calculate SNP Index

#'
#' This function is used to sparate (filter) the two bulks, and calculate their (SNP index) as follows:
#' \deqn{SNPindex=as.numeric(AD_alt)/(as.numeric(AD_ref) + as.numeric(AD_alt)}. 
#' Bulks are finally joined together in a single data frame.
#'
#' @param vcf.df Data frame of the vcf file
#' @param WTbulk Wild-Type pool
#' @param Mbulk Mutant pool 
#'
#' @return A dataframe with the joint bulks
#'
#' @export
#' @examples
#' calc_SNPindex()

calc_SNPindex <- function(vcf.df, WTbulk, Mbulk) {

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