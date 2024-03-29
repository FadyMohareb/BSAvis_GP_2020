% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_SNPratio.R
\name{calc_SNPratio}
\alias{calc_SNPratio}
\title{Calculate SNP-ratio}
\usage{
calc_SNPratio(vcf.df, wtBulk, mBulk, variants = "SNP")
}
\arguments{
\item{vcf.df}{Data frame of the vcf file}

\item{wtBulk}{Wild-Type pool}

\item{mBulk}{Mutant pool}

\item{variants}{variants to be considered. Default is "SNP" (allowed: "SNP" or "all")}
}
\value{
Data frame containing filtered variant and SNP-ratio information.
}
\description{
This function allows to calculate the SNP-ratio between two bulks (M/WT), as follows:

\deqn{(AD_ref.WT/(AD_ref.WT + AD_alt.WT)) - (AD_ref.M/(AD_ref.M + AD_alt.M))}

Calculated SNP-ratio values get added in a new column of the data frame (returned by the readBSA_vcf() function).
Note that the user can select specific variants to consider, by setting the "Variants" parameter to "SNP" (default) or "all" (respectively referring to SNPs or InDels+SNPs).
}
\details{
The data frame returned by readBSA_vcf() is filtered by bulk ID to create two separate data frames: 
one specific to the wild-type bulk variants information and other one specific to the mutant bulk variants information.
In the case the type of variants is set to "SNP", variants corresponding to InDels are discarded by removing the rows in the data frame containing more than three characters in the GT_alleles column (e.g.,"T/AT"corresponds to an insertion and contains 4 characters). 
This also applies for deletions, where one of the three characters is a "*". Then, both data frames (wild-type bulk and mutant bulk) are joined by rows having the same value in the chromosome and position columns. 
The SNP-ratio is then calculated for each of the variants and the values are added in a new column of the joint data frame, which will be returned by the function.
}
\examples{
## Calculate SNP-index for both bulks (only SNPs will be considered)
vcf_df_SNPratio <- calc_SNPratio(vcf.df=vcf_list$df, 
                                 wtBulk="pool_S3781_minus", 
                                 mBulk="pool_S3781_plus") 
## Calculate SNP-index considering both InDels and SNPs
vcf_df_SNPratio <- calc_SNPratio(vcf.df=vcf_list$df, 
                                 wtBulk="pool_S3781_minus", 
                                 mBulk="pool_S3781_plus", 
                                 variants="all")
}
\references{
Wachsman et al., 2017
}
