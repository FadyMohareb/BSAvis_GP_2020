![BSAvis logo](BSAvis_small_logo.png)<br>
BSAvis is a flexible, user-friendly package built for Bulk Segregant Analyses (BSA), capable of generating publication-quality plots to visualise and identify possible causal regions related to traits of interest in bulks expressing contrasting phenotypes.

## Table of Contents
> - [Prerequisites](#Prerequisites)<br>
> - [Documentation](#Documentation)<br>
> - [Testing Files](#Testing-files)<br>
> - [BSAvis Dependencies](#BSAvis-Dependencies)<br>
> - [Installing BSAvis](#Installing-BSAvis)<br>
> - [Reading the VCF file](#Reading-the-VCF-file)<br><br>
> - [BSAvis Package - Practical Examples](#BSAvis-Package-\--Practical-Examples)<br>
>   - [Combined BSA and Plotting](#Combined-BSA-and-Plotting)<br>
>       - [SNP-index Method](#snp-index-method)<br>
>       - [delta(SNP-index) Method](#deltasnp-index-method)<br>
>       - [SNP-ratio Method](#snp-ratio-method)<br>
>   - [Stepwise BSA and Plotting](#Stepwise-BSA-and-Plotting)<br>
>       - [SNP-index/∆(SNP-index) Method](#snp-indexsnp-index-method)<br>
>       - [SNP-ratio Method](#snp-ratio-method-1)<br><br>
> - [Interactive BSAvis Dashboard](#Interactive-BSAvis-Dashboard)
>   - [Installing Shiny Libraries](#installing-shiny-libraries)<br>
>   - [Running BSAvis R-Shiny Application](#running-bsavis-r-shiny-application)<br><br>
> - [Authors](#Authors)


# Prerequisites
BSAvis is compatible for being run on either MacOS or Windows operating system.<br>
R and RStudio integrated development environment (IDE) (version >=3.6.1) are required to run BSAvis. <br>

- [RStudio](https://www.rstudio.com/products/rstudio/download/ "RStudio") (free version) <br>
- [R](https://cran.r-project.org "R") 

BSAvis package requires merged Variant Calling Format (`VCF`) files as input files, generated using `GATK4`.<br> For best results from the BSAvis package, the following sample scripts are recommended:
- [Alignment Steps](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/QC_Alignment_VC/alignment_variantCalling/steps/alignment_steps.txt "Alignment Steps")
- [Variant Calling Steps](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/QC_Alignment_VC/alignment_variantCalling/steps/variantCalling_steps.txt "Variant Calling Steps")

Be aware that the `AD` (Allelic Depth) column needs to be included inside the `VCF` file (obtained using `GATK4`). Previous versions to this are, therefore, not recommended. <br> Joint genotyping is also required to obtain a single `VCF` file that includes variants from both bulks. 

**This tool was tested and run both on Windows and MacOS.**

# Documentation
A brief step-by-step process on how to run BSAvis is included on this page, using testing files. However, please refer to the [User Manual](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/Documentation/user_manual_BSAvis.pdf "BSAvis User Manual") and [Technical Documentation](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/Documentation/technical_documentation_BSAvis.pdf "Technical Documentation") for a better understanding of the implemented functions.

# Test Files
The following testing files have been provided:
* `test.RData`
* `test.vcf` 

Note that `test.RData` can be downloaded directly by clicking on the "download" button ([on this page](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/test.RData)). The actual `VCF` file can be obtained only after downloading and unzipping the full repository (click on the green "Code" button at the top-right corner of this repository and select "Download ZIP"). <br> Refer to the [User Manual](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/Documentation/user_manual_BSAvis.pdf "BSAvis User Manual") for guidance.

# BSAvis Dependencies
Before downloading the BSAvis package and to avoid encountering issues, it is strongly recommended to ensure having the following libraries installed and loaded on RStudio: `vcfR`, `ggplot2`, `dplyr`, `tidyr`and `devtools`.

```R
# Install packages 
install.packages(c("vcfR", "ggplot2", "dplyr", "tidyr", "devtools"))

# Load packages
library(vcfR)
library(ggplot2)
library(dplyr)
library(tidyr)
```

# Installing BSAvis
To install the BSAvis package, run the following commands:

```R
#Install BSAvis from GitHub
devtools::install_github("FadyMohareb/BSAvis_GP_2020/BSAvis")

#Load BSAvis
library(BSAvis)
```
Please note that **help pages** for every implemented function can be inspected on RStudio at any time:

```R
#Common
?readBSA_vcf 
#BSAvis Package Functions
?SNPindex_plot
?deltaSNPindex_plot  
?SNPratio_plot 
?calc_SNPindex 
?calc_SNPratio 
?calc_deltaSNPindex 
?extract_chrIDs 
?filter_SNPindex 
?filter_SNPratio 
?plot_SNPindex 
?plot_SNPratio 
?plot_deltaSNPindex 
#BSAvis R-Shiny App
?BSAvis_shiny
#NOT to be run manually
?shiny_SNPindex 
?shiny_SNPratio 
?shiny_deltaSNPindex 
?slidingWindow
?shinyPlot_SNPindex 
?shinyPlot_deltaSNPindex 
?shinyPlot_SNPratio 
```

# Reading the VCF file
After successfully loading the BSAvis package, the merged `VCF` file needs to be read and loaded inside the working environment. Be aware that **this step is common for every implemented method** and is needed to run the rest of the package functions. 

Proceed running the following command:

```R
#Read VCF file
vcf_list <- readBSA_vcf("test.vcf")
```

**Warning –** this step might take time to complete (approximately 15-20 minutes). It is strongly recommended to leave RStudio open and running, to properly process the data. Move to the next step only after the red button in the top right corner of the console disappears.

Alternatively, to familiarize yourself with the package and skip the `readBSA_vcf()`function, you can directly load the `vcf_list` object using the provided `.RData` file (`test.RData`), either by double-clicking on the file or by running the following command:

```R
#Load workspace
load("test.RData")
```

# BSAvis Package - Practical Examples
# Combined BSA and Plotting
Wrapper functions included in this section will apply the chosen BSA method and return plots in `TIFF` format. 

A step-by-step approach to BSA and plotting is also included (see next section).

### SNP-index Method
```R
#SNP-index Method
SNPindex_plot(vcf.list=vcf_list, 
wtBulk="pool_S3781_minus",
mBulk="pool_S3781_plus",
variants="SNP",
min.SNPindex=0.3,
max.SNPindex=0.9,
min.DP=50,
max.DP=200,
min.GQ=99,
chrID="SL4.0ch03",
chr=3,
windowSize=1000000,
windowStep=10000, 
filename="plot_SNPindex_ch", 
path="/currentWorkingDirectory", 
dpi=1200, # if set, the plot gets saved 
width=7.5,
height=5, 
units="in")
```
### delta(SNP-index) Method
```R
#delta(SNP-index) Method
deltaSNPindex_plot(vcf.list=vcf_list, 
wtBulk="pool_S3781_minus",
mBulk="pool_S3781_plus",
variants="SNP",
min.SNPindex=0.3,
max.SNPindex=0.9,
min.DP=50,
max.DP=200,
min.GQ=99,
chrID="SL4.0ch03",
chr=3,
windowSize=1000000,
windowStep=10000, 
filename="plot_deltaSNPindex_ch", 
path="/currentWorkingDirectory", 
dpi=1200, # if set, the plot gets saved 
width=7.5,
height=5, 
units="in")
```
### SNP-ratio Method
```R
#SNP-ratio Method
SNPratio_plot(vcf.list=vcf_list, 
wtBulk="pool_S3781_minus",
mBulk="pool_S3781_plus",
variants="SNP", 
min.SNPratio=0.1, 
min.DP=50, 
max.DP=200, 
chrID="SL4.0ch03", 
chr=3,
degree=2,
span=0.03, 
filename="plot_SNPratio_ch", 
path="/currentWorkingDirectory",
dpi=1200, # if set, the plot gets saved 
width=7.5,
height=5,
units="in")
```
### Important Notes
For all functions listed above, please note that:

- The minimum required parameters are: `vcf.list`, `wtBulk`, `mBulk`, `chrID` and `chr`
- `filename`, `path`, `width`, `height` and `units` are all part of the plot-saving functionality and are directly linked with the `dpi` parameter. Without setting `dpi`, all the previously mentioned ones will be ignored
- Some parameters are already set to default but can be customized
- To properly implement the function, please refer to the [Technical Documentation](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/Documentation/technical_documentation_BSAvis.pdf "Technical Documentation") or inspect the [help pages](#Installing-BSAvis) of each function on RStudio 

# Stepwise BSA and Plotting
Functions included in this section will apply the chosen BSA method and return the plots to the user in a multistep process, conversely to the previous section which describes the process for a simple combined functionality approach to BSA and plotting.<br>

All of the included examples refer have been applied to the **testing dataset**.<br>

**Important notes** 
- be sure to have run the `readBSA_vcf()`function before moving on
- some parameters are already set to default but can be customized
- to properly implement the functions, please refer to the [Technical Documentation](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/Documentation/technical_documentation_BSAvis.pdf "Technical Documentation") or directly inspect the help pages on RStudio.

## SNP-index/∆(SNP-index) Method
1.	Calculate SNP-indices for both bulks using the `calc_SNPindex()` function
Example:
```R
vcf_df_SNPindex <- calc_SNPindex(vcf.df=vcf_list$df,
                    wtBulk="pool_minus", 
                    mBulk="pool_plus",
                    variants="SNP")
 ```


2.	Filter variants using the `filter_SNPindex()` function.

Example:
```R
vcf_df_SNPindex_filt <- filter_SNPindex(
                        vcf.df.SNPindex=vcf_df_SNPindex, 
                        min.SNPindex=0.3, 
                        max.SNPindex=0.9, 
                        min.DP=50, 
                        max.DP=200, 
                        min.GQ=99)
```

3.	Extract chromosome IDs using the `extract_chrIDs()` function.

Example:
```R
chromList <- extract_chrIDs(meta=vcf_list$meta)
```
4.	Calculate the sliding windows based on the chromosome length of the specified  
chromosome using the `slidingWindow()` function.

Example:
```R
SNPindex_windows <- slidingWindow(meta=vcf_list$meta, 
                    chrList=chromList, 
                    chrID="SL4.0ch03", 
                    windowSize=1000000, 
                    windowStep=10000, 
                    vcf.df.SNPindex.filt=vcf_df_SNPindex_filt)
```
5.	Plot SNP-index across the positions of a given chromosome using the  
`plot_SNPindex()` function.

Example:
```R
plot_SNPindex(SNPindex.windows=SNPindex_windows,
chr=3,
filename="plot_SNPindex_ch",
path="currentWorkingDirectory",
dpi=1200,
width=7.5,
height=5,
units="in") 
```

6.	Calculate ∆(SNP-index) using the `calc_deltaSNPindex()` function.

Example:
```R
deltaSNPindex_windows <- calc_deltaSNPindex(SNPindex.windows=SNPindex_windows)
```

7.	Plot ∆(SNP-index) across the positions of a given chromosome using the `plot_deltaSNPindex()` function.

Example:
```R
plot_deltaSNPindex(deltaSNPindex.windows=deltaSNPindex_windows,
chr=3,
filename="plot_deltaSNPindex_ch",
path="currentWorkingDirectory", 
dpi=1200, 
width=7.5,
height=5,
units="in")
```
 
## SNP-ratio Method
1.	Calculate SNP-ratios for both bulks using the `calc_SNPratio()` function

Example:
```R
vcf_df_SNPratio <- calc_SNPratio(vcf.df=vcf_list$df, 
                    wtBulk="pool_S3781_minus", 
                    mBulk="pool_S3781_plus", 
                    variants="SNP")
```

2.	Filter variants using the `filter_SNPratio()` function

Example:
```R
vcf_df_SNPratio_filt <- filter_SNPratio(
                        vcf.df.SNPratio=vcf_df_SNPratio, 
                        min.SNPratio=0.1, 
                        min.DP=50, 
                        max.DP=200)
```
3.	Extract chromosome IDs using the `extract_chrIDs()`function

Example:
```R
chromList <- extract_chrIDs(vcf_list$meta)
```

4.	Plot SNP-ratio across the positions of a given chromosome using the `plot_SNPratio()` function

Example:
```R
plot_SNPratio(vcf.df.SNPratio.filt=vcf_df_SNPratio_filt, 
chrList=chromList, 
chrID="SL4.0ch03", 
chr=3, 
min.SNPratio=0.1, 
degree=2, 
span=0.3, 
filename="plot_SNPratio_ch", 
path="currentWorkingDirectory", 
dpi=1200, 
width=7.5, 
height=5, 
units="in")
```

# Interactive BSAvis Dashboard
This section describes how to run the interactive version of BSAvis, which allows the user to perform interactive BSA analysis using a user-friendly R-Shiny application.

## Installing Shiny Libraries
To begin, the user is required to manually install and load the required libraries on RStudio, following two steps:

1.	Install the required libraries using the following command: 

```R
install.packages(c("shiny", "shinycssloaders", "shinyalert"))
```

Be patient and proceed only after the installation is completed.

2.	Load the libraries running the following commands: 
```R
library(shiny)
library(shinycssloaders)
library(shinyalert)
```

## Running BSAvis R-Shiny Application
BSAvis interactive tool can be run by calling the `BSAvis_shiny()` function. 

Be sure to have run the `readBSA_vcf()` function, as recommended at the beginning, since its output (`vcf_list`) is needed for the `BSAvis_shiny()` function:
```R
BSAvis_shiny(vcf_list)
```

**Important Notes**
- All generated plots can be saved at any time, by clicking on the `"Save Plot"` button, found under each plotting panel
- Plots can be zoomed in by dragging, releasing, and double-clicking on the selected area
- Please refer to the [User Manual](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/Documentation/user_manual_BSAvis.pdf "BSAvis User Manual") for a better understanding of the BSAvis R-Shiny application, together with its functionalities.

# Authors
BSAvis (Version 1.0) was developed by:<br>

_Elisabetta Galatola, Rebecca Guy, Weiyi Huang, Sweta Jajoriya, Claudia Rey-Carrascosa_<br><br>
_Applied Bioinformatics MSc_<br>
_Cranfield University - Cranfield, Bedford, UK_<br>
_Academic Year: 2020-2021_