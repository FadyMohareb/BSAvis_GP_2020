# Alignment
Each of the Illumina Pair-End reads were aligned to the reference tomato genome (*Solanum lycopersicum, version 4.0* - downloaded from the Solanaceae Genomics Network: https://solgenomics.net), using BWA aligner (version 0.7.7). 

The resulting **SAM** file was sorted by coordinates and converted to **BAM** format using Picard tools (*SortSam* function). PCR duplicates were also marked using Picard tools (*MarkDuplicates* function).

For best results from the BSAvis package, the following sample script is recommended: [alignment_steps.txt](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/QC_Alignment_VC/alignment_variantCalling/steps/alignment_steps.txt)

# Variant Calling
Variant calling was performed for each sample using *HaplotypeCaller* function from the Genome Analysis Toolkit (GATK4, version 4.1.9.0), in GVCF mode. Single-sample GVCFs generated from both datasets (i.e., two parents and two pools in each dataset) were then imported into GenomicsDB and merged using the *GenomicsDBImport* function. Finally, GATK4 *GenotypeGVCFs* tool was used to perform joint genotyping on the dataset samples, previously loaded into the GenomicsDB workspace. As a result, a single VCF file per dataset was generated, containing variant information of all four samples.

For best results from the BSAvis package, the following sample script is recommended: [alignment_steps.txt](https://github.com/FadyMohareb/BSAvis_GP_2020/blob/main/QC_Alignment_VC/alignment_variantCalling/steps/variantCalling_steps.txt) 