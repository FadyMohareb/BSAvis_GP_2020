#Extract Chromosome IDs

#' @title Extract Chromosome IDs
#' @description This function allows to extract the chromosome IDs, based on the meta information contained in the vcf file.
#'
#' @param meta meta information 
#'
#' @return Chromosome list 
#'
#' @export extract_chrIDs
#' @examples
#' chromList <- extract_chrIDs(meta=vcf_list$meta)


extract_chrIDs <- function(meta) {
  
  #Extract those meta lines containing the length of the chromosomes
  lengthLines <- meta[grep("<ID.*length=", meta)]
  #Remove characters after chromosome ID (by replacing with "")
  chrList <- sub(",length.*", "", lengthLines)
  #Remove characters before chromosome ID (by replacing with "")
  chrList <- sub(".*ID=", "", chrList) #Contains all chromosomes IDs
  
  return(chrList)
}
