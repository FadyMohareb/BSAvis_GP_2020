# BSAvis Package
Please refer to the [User Manual]() and [Technical Documentatio]() for a better understanding of the package.

```R
#Install BSAvis package
devtools::install_github("FadyMohareb/BSAvis_GP_2020/BSAvis")
library(BSAvis)

#Inspect functions
#Common
?readBSA_vcf 
#BSAvis Package Functions
#Combined BSA and Plotting (wrapper functions)
?SNPindex_plot
?deltaSNPindex_plot  
?SNPratio_plot 
#Stepwise BSA and Plotting
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