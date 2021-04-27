#Install bsaR package
remotes::install_github("xuzhougeng/bsaR")

#Load additional packages
library(vcfR)
library(bsaR)

#Store vcf file name
vcf.file <- "dataset2.vcf"

#Choose which chromosomes will be included in the analysis 
Chroms <- paste0(rep("SL4.0ch0", 9), 1:9)
Chroms <- append(Chroms, paste0(rep("SL4.0ch", 3), 10:12))

#Create the BSA project from the vcf file
bsa <- CreateBsaFromVcf(vcf.file = vcf.file)

#Filter the BSA 
bsa <- FilterMultiVariant(bsa)
bsa <- FilterNaVariant(bsa)

bsa <- CalcDepth(bsa)
bsa <- FilterByDepth(bsa, min.depth = 40)

#Calculate the frequency
bsa <- CalcAltFreq(bsa)
bsa <- CalcFreqByWindow(bsa, window.size = 700000)

#Visualise plot
pdf("plots_bsaR_deltaSNP_ds2.pdf", title="Delta(SNP-index) S3785 pools")
for (i in Chroms){
  plotWindowFreq(bsa, i, "pool_S3785_plus", "pool_S3785_minus",
                 window.size = 700000)
  mtext(paste("Delta(SNP-index) S3785 pools - Chromosome", i), outer=TRUE,  cex=1, line=-3)
}
dev.off()
