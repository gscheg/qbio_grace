#Final Project code

#mutation data vs. transcriptomics

if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
#BiocManager::install("maftools")
library(maftools)
library(SummarizedExperiment)

#load in clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

#subtypes <- TCGAquery_subtype(tumor = "BRCA")

#MAF
mutation <- GDCquery_Maf(tumor = "BRCA", save.csv=TRUE, pipeline="varscan2")
colnames(clinic)[1] <- "Tumor_Sample_Barcode"
my_maf <- data.table::fread("GDCdata/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.csv")
#reads in as maf
maf <- read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

TP53mut <- subsetMaf(maf, gene = "TP53")
PIK3CAmut <- subsetMaf(maf, gene = "PIK3CA")
TTNmut <- subsetMaf(maf, gene = "TTN")

#maf@clinical.data$Tumor_Sample_Barcode
TP53mutbarcodes <- TP53mut@clinical.data$Tumor_Sample_Barcode #patients with mutated TP53
PIK3CAmut@clinical.data$Tumor_Sample_Barcode
TTNmut@clinical.data$Tumor_Sample_Barcode

#barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
                     #"TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
                     #"TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
#load in HTSeq data
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts")
GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

#Get counts for specific genes
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient #shorter barcodes to match those in MAF

#Access TP53 counts
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TP53_mask]
mut_TP53_counts <- htseq_counts[TP53_ENSG_name, TP53mutbarcodes]
nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, !(colnames(htseq_counts) %in% TP53mutbarcodes)]
boxplot(mut_TP53_counts, nonmut_TP53_counts)

