#Install QTLseqr
install.packages("devtools")
devtools::install_github("bmansfeld/QTLseqr")
#Intsall additional needed package
BiocManager::install("genefilter")

#load the package
library("QTLseqr")

#Choose which chromosomes will be included in the analysis 
Chroms <- paste0(rep("SL4.0ch0", 9), 1:9)
Chroms <- append(Chroms, paste0(rep("SL4.0ch", 3), 10:12))

#Import SNP data from file
df <- importFromGATK(
  file = "dataset1_pools_QTLseqr.tsv",
  highBulk = "pool_S3781_plus",
  lowBulk = "pool_S3781_minus",
  chromList = Chroms
)

#Filter SNPs based on some criteria
df_filt <-  filterSNPs(
  SNPset = df,
  refAlleleFreq = 0.20,
  minTotalDepth = 50,
  maxTotalDepth = 400,
  minSampleDepth = 40,
  minGQ = 99
)

#Run G' analysis
df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP"
)

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  popStruc = "F2",
  bulkSize = c(29, 39),
  replications = 10000,
  intervals = c(95, 99)
)

pdf("C:\\Users\\s341924\\OneDrive - Cranfield University\\GP\\Plots\\plots_QTLseqR_ds1.pdf", title="S3781 pools")
for (Chrom in Chroms){
  #Plot Gprime
  Gprime <- plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01, subset=Chrom)
  print(Gprime)
}
for (Chrom in Chroms){
  #Plot delta SNP-index
  SNPindex <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE, subset=Chrom)
  print(SNPindex)
}
for (Chrom in Chroms){
  #Plot nSNPs
  nSNPs <- plotQTLStats(SNPset = df_filt, var = "nSNPs", subset=Chrom)
  print(nSNPs)
}
for (Chrom in Chroms){
  #Plot negLog10Pval
  negLog10Pval <- plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, subset=Chrom)
  print(negLog10Pval)
}
dev.off()


#Plot all chromosomes in a single graph
#Plot Gprime
png("plot_QTLseqr_Gprime_ds1.png", width=1300, height=500)
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
dev.off()

#Plot delta SNP-index
png("plot_QTLseqr_deltaSNP_ds1.png", width=1300, height=500)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
dev.off()

#Plot nSNPs
png("plot_QTLseqr_nSNPs_ds1.png", width=1300, height=500)
plotQTLStats(SNPset = df_filt, var = "nSNPs")
dev.off()

#Plot negLog10Pval
png("plot_QTLseqr_negLog10Pval_ds1.png", width=1300, height=500)
plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01)
dev.off()

