#Extract Chromosome Lenght

#'
#' This function allows to extract the chromosome length of a specific chromosome, using the meta information contained in the vcf file.
#'this function extract chromosome length
#'
#' @param chromID chromosome ID 
#' @param vcf
#'
#' @export
#' @examples
#' extract_chromLength()

extract_chromLength <- function(chromID, meta){
  
  #Extract meta line containing length info on the specific chromosome
  lengthLine <- meta[grep(chromID, meta)]
  #Remove (by replacing with "") previous characters to length number
  length <- sub(".*length=", "", lengthLine)
  #Set chromosome size (and remove ">" character which appears after length number)
  chromLength <- as.numeric(sub(">", "", length))
  
  return(chromLength)
}
