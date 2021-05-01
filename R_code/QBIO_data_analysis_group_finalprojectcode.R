#Final Project code

#boxplots of transcriptomics vs. mutation data

if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
#BiocManager::install("maftools")
library(maftools)
library(SummarizedExperiment)

#load in clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
#GDCdownload(clin_query) #only need once
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")

#load in subtypes info
subtypes <- TCGAquery_subtype(tumor = "BRCA")

#barcodes for each subtype
LumA_mask <- subtypes[subtypes$"BRCA_Subtype_PAM50"=="LumA", ]
LumAbarcodes <- LumA_mask$patient
LumB_mask <- subtypes[subtypes$"BRCA_Subtype_PAM50"=="LumB", ]
LumBbarcodes <- LumB_mask$patient
Basal_mask <- subtypes[subtypes$"BRCA_Subtype_PAM50"=="Basal", ]
Basalbarcodes <- Basal_mask$patient
Her2_mask <- subtypes[subtypes$"BRCA_Subtype_PAM50"=="Her2", ]
Her2barcodes <- Her2_mask$patient
Normal_mask <- subtypes[subtypes$"BRCA_Subtype_PAM50"=="Normal", ]
Normalbarcodes <- Normal_mask$patient

#load in MAF
mutation <- GDCquery_Maf(tumor = "BRCA", save.csv=TRUE, pipeline="varscan2")
colnames(clinic)[1] <- "Tumor_Sample_Barcode"
my_maf <- data.table::fread("GDCdata/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.csv")
#reads in as maf
maf <- read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

#subset of MAF data for each gene of interest
TP53mut <- subsetMaf(maf, gene = "TP53")
PIK3CAmut <- subsetMaf(maf, gene = "PIK3CA")

#barcodes of patients with mutated gene for both genes of interest
TP53mutbarcodes <- TP53mut@clinical.data$Tumor_Sample_Barcode #patients with mutated TP53
PIK3CAmutbarcodes <- PIK3CAmut@clinical.data$Tumor_Sample_Barcode

#load in HTSeq data
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts", barcode = barcodes_rnaseq)
#GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient #shorter barcodes to match those in MAF

#barcodes of patients in each subtype in the HTSeq data
LumA <- colnames(htseq_counts) %in% LumAbarcodes
LumB <- colnames(htseq_counts) %in% LumBbarcodes
Basal <- colnames(htseq_counts) %in% Basalbarcodes
Her2 <- colnames(htseq_counts) %in% Her2barcodes
Normal <- colnames(htseq_counts) %in% Normalbarcodes

#access TP53 counts
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TP53_mask]
#access mutated TP53 barcodes in HTSeq data
TP53mutated <- colnames(htseq_counts) %in% TP53mutbarcodes

#access PIK3CA counts
PIK3CA_mask <- rowData(sum_exp)$external_gene_name == "PIK3CA"
PIK3CA_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[PIK3CA_mask]
#access mutated PIK3CA barcodes in HTSeq data
PIK3CAmutated <- colnames(htseq_counts) %in% PIK3CAmutbarcodes

#TP53
#get log counts for each category: each combination of TP53 mutational status 
#(mutated or non-mutated) and the five subtypes
LumA_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, LumA & TP53mutated]
if(length(LumA_mut_TP53_counts) != 0){ #if the counts is not 0, can apply log10
  LumA_mut_TP53_counts <- log10(LumA_mut_TP53_counts)
}
LumA_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, LumA & !TP53mutated]
if(length(LumA_nonmut_TP53_counts) != 0){
  LumA_nonmut_TP53_counts <- log10(LumA_nonmut_TP53_counts)
}
LumB_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, LumB & TP53mutated]
if(length(LumB_mut_TP53_counts) != 0){
  LumB_mut_TP53_counts <- log10(LumB_mut_TP53_counts)
}
LumB_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, LumB & !TP53mutated]
if(length(LumB_nonmut_TP53_counts) != 0){
  LumB_nonmut_TP53_counts <- log10(LumB_nonmut_TP53_counts)
}
Basal_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, Basal & TP53mutated]
if(length(Basal_mut_TP53_counts) != 0){
  Basal_mut_TP53_counts <- log10(Basal_mut_TP53_counts)
}
Basal_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, Basal & !TP53mutated]
if(length(Basal_nonmut_TP53_counts) != 0){
  Basal_nonmut_TP53_counts <- log10(Basal_nonmut_TP53_counts)
}
Her2_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, Her2 & TP53mutated]
if(length(Her2_mut_TP53_counts) != 0){
  Her2_mut_TP53_counts <- log10(Her2_mut_TP53_counts)
}

Her2_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, Her2 & !TP53mutated]
if(length(Her2_nonmut_TP53_counts) != 0){
  Her2_nonmut_TP53_counts <- log10(Her2_nonmut_TP53_counts)
}

Normal_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, Normal & TP53mutated]
if(length(Normal_mut_TP53_counts) != 0){
  Normal_mut_TP53_counts <- log10(Normal_mut_TP53_counts)
}

Normal_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, Normal & !TP53mutated]
if(length(Normal_nonmut_TP53_counts) != 0){
  Normal_nonmut_TP53_counts <- log10(Normal_nonmut_TP53_counts)
}

pdf("TP53_RNAmut_boxplot.pdf")
par(mar = c(9,4,4,2) + 0.1) #fits labels in boxplot display
boxplot(LumA_mut_TP53_counts, LumA_nonmut_TP53_counts, LumB_mut_TP53_counts, LumB_nonmut_TP53_counts, Basal_mut_TP53_counts, 
        Basal_nonmut_TP53_counts, Her2_mut_TP53_counts, Her2_nonmut_TP53_counts, Normal_mut_TP53_counts, Normal_nonmut_TP53_counts,
        names = c("LumA mut", "LumA non-mut", "LumB mut", "LumB non-mut", "Basal mut", "Basal non-mut",
                  "Her2 mut", "Her2 non-mut", "Normal mut", "Normal non-mut"), 
        ylab = "HTSeq mRNA log counts for TP53 gene", las = 2)
dev.off()

#PIK3CA
#get log counts for each category: each combination of PIK3CA mutational status 
#(mutated or non-mutated) and the five subtypes
LumA_mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, LumA & PIK3CAmutated]
if(length(LumA_mut_PIK3CA_counts) != 0){
  LumA_mut_PIK3CA_counts <- log10(LumA_mut_PIK3CA_counts)
}
LumA_nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, LumA & !PIK3CAmutated]
if(length(LumA_nonmut_PIK3CA_counts) != 0){
  LumA_nonmut_PIK3CA_counts <- log10(LumA_nonmut_PIK3CA_counts)
}
LumB_mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, LumB & PIK3CAmutated]
if(length(LumB_mut_PIK3CA_counts) != 0){
  LumB_mut_PIK3CA_counts <- log10(LumB_mut_PIK3CA_counts)
}
LumB_nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, LumB & !PIK3CAmutated]
if(length(LumB_nonmut_PIK3CA_counts) != 0){
  LumB_nonmut_PIK3CA_counts <- log10(LumB_nonmut_PIK3CA_counts)
}
Basal_mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, Basal & PIK3CAmutated]
if(length(Basal_mut_PIK3CA_counts) != 0){
  Basal_mut_PIK3CA_counts <- log10(Basal_mut_PIK3CA_counts)
}
Basal_nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, Basal & !PIK3CAmutated]
if(length(Basal_nonmut_PIK3CA_counts) != 0){
  Basal_nonmut_PIK3CA_counts <- log10(Basal_nonmut_PIK3CA_counts)
}
Her2_mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, Her2 & PIK3CAmutated]
if(length(Her2_mut_PIK3CA_counts) != 0){
  Her2_mut_PIK3CA_counts <- log10(Her2_mut_PIK3CA_counts)
}
Her2_nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, Her2 & !PIK3CAmutated]
if(length(Her2_nonmut_PIK3CA_counts) != 0){
  Her2_nonmut_PIK3CA_counts <- log10(Her2_nonmut_PIK3CA_counts)
}

Normal_mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, Normal & PIK3CAmutated]
if(length(Normal_mut_PIK3CA_counts) != 0){
  Normal_mut_PIK3CA_counts <- log10(Normal_mut_PIK3CA_counts)
}
Normal_nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, Normal & !PIK3CAmutated]
if(length(Normal_nonmut_PIK3CA_counts) != 0){
  Normal_nonmut_PIK3CA_counts <- log10(Normal_nonmut_PIK3CA_counts)
}

pdf("PIK3CA_RNAmut_boxplot.pdf")
par(mar = c(9,4,4,2) + 0.1)
boxplot(LumA_mut_PIK3CA_counts, LumA_nonmut_PIK3CA_counts, LumB_mut_PIK3CA_counts, LumB_nonmut_PIK3CA_counts, Basal_mut_PIK3CA_counts,
        Basal_nonmut_PIK3CA_counts, Her2_mut_PIK3CA_counts, Her2_nonmut_PIK3CA_counts, Normal_mut_PIK3CA_counts, Normal_nonmut_PIK3CA_counts,
        names = c("LumA mut", "LumA non-mut", "LumB mut", "LumB non-mut", "Basal mut", "Basal non-mut",
                  "Her2 mut", "Her2 non-mut", "Normal mut", "Normal non-mut"), 
        ylab = "HTSeq mRNA log counts for PIK3CA gene", las = 2)
dev.off()
