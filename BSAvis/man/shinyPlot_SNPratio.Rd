% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyPlot_SNPratio.R
\name{shinyPlot_SNPratio}
\alias{shinyPlot_SNPratio}
\title{Shiny: Plot SNP-ratio}
\usage{
shinyPlot_SNPratio(
  vcf.df.SNPratio.filt,
  chrList,
  chrID,
  chr,
  min.SNPratio,
  degree = 2,
  span = 0.07,
  ranges
)
}
\arguments{
\item{vcf.df.SNPratio.filt}{filtered data frame (containing both bulks)}

\item{chrList}{list of chromosome IDs}

\item{chrID}{chromosome ID of interest}

\item{chr}{chromosome number to print on plot}

\item{min.SNPratio}{min SNP ratio threshold}

\item{degree}{LOESS smoothing degree (default=2)}

\item{span}{LOESS smoothing span (default=0.07)}

\item{ranges}{axes ranges (x,y)}
}
\description{
Please note that this function is not to be run manually.
}