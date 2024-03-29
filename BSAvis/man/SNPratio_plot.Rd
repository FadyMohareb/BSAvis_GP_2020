% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNPratio_plot.R
\name{SNPratio_plot}
\alias{SNPratio_plot}
\title{SNP-ratio Wrapper Function}
\usage{
SNPratio_plot(
  vcf.list,
  wtBulk,
  mBulk,
  variants = "SNP",
  min.SNPratio = 0.1,
  min.DP = 50,
  max.DP = 200,
  chrID,
  chr,
  degree = 2,
  span = 0.07,
  filename = paste0("plot_SNPratio_ch", chr),
  path = getwd(),
  dpi,
  width = 7.5,
  height = 5,
  units = "in"
)
}
\arguments{
\item{vcf.list}{object containing meta information and vcf data frame}

\item{wtBulk}{Wild-Type pool}

\item{mBulk}{Mutant pool}

\item{variants}{variants to be considered. Default is "SNP" (allowed: "SNP" or "all")}

\item{min.SNPratio}{min value allowed for the SNP index (default=0.3)}

\item{min.DP}{min value allowed for the read depth (default=50)}

\item{max.DP}{max value allowed for the read depth (default=200)}

\item{chrID}{chromosome ID of interest}

\item{chr}{chromosome name printed on the plot}

\item{degree}{LOESS smoothing degree (default=2)}

\item{span}{LOESS smoothing span (default=0.07)}

\item{filename}{file name under which the file will be saved (default="plot_SNPratio_chX")}

\item{path}{path where the file will be saved (default=current working directory)}

\item{dpi}{resolution value. If no value is given, plots will be generated but not saved}

\item{width}{width value (default=7.5)}

\item{height}{height value (default=5)}

\item{units}{size units (default="in")}
}
\description{
This wrapper function is used to fully run the SNP-ratio method, by calling all the funcitons involved in plotting the SNP-ratio values across positions of a specific chromosome.

\deqn{SNPindex=AD_alt/(AD_ref + AD_alt)}
}
\details{
Wrapper function that sequentially calls the required functions involved in generating the SNP-ratio plot:
calc_SNPratio(), filter_SNPratio(), extract_chrIDs() and plot_SNPratio().
The resulting plot will show SNP-ratio values (for both bulks) across chromosome positions.
}
\examples{
## Use default values WITHOUT saving the plot
SNPratio_plot(vcf.list=vcf_list, 
         wtBulk="pool_S3781_minus", 
         mBulk="pool_S3781_plus", 
         chrID="SL4.0ch03", 
         chr=3)
## OR use default values AND save the plot
SNPratio_plot(vcf.list=vcf_list, 
         wtBulk="pool_S3781_minus", 
         mBulk="pool_S3781_plus", 
         chrID="SL4.0ch03", 
         chr=3, 
         dpi=1200)
## OR customise default parameters
SNPratio_plot(vcf.list=vcf_list, 
         wtBulk="pool_S3781_minus", 
         mBulk="pool_S3781_plus", 
         variants="all",
         min.SNPratio=0.3, 
         min.DP=60, 
         max.DP=250,
         chrID="SL4.0ch03", 
         chr=3, 
         degree=1, 
         span=0.05, 
         filename="deltaSNPindex_chrom03", 
         path="Document/Plots", 
         dpi=1200, 
         width=20, 
         height=12, 
         units="cm")
}
