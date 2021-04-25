#Final Project code

#mutation data vs. transcriptomics

if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
#BiocManager::install("maftools")
library(maftools)
library(SummarizedExperiment)

#load in clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

subtypes <- TCGAquery_subtype(tumor = "BRCA")
age_subs = subtypes$age_at_initial_pathologic_diagnosis
#also add the age category to the subtypes information
subtypes$age_category = ifelse(age_subs < 40, "Young", ifelse(age_subs >= 60, "Old", "Mid"))

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

#MAF
#mutation <- GDCquery_Maf(tumor = "BRCA", save.csv=TRUE, pipeline="varscan2")
colnames(clinic)[1] <- "Tumor_Sample_Barcode"
my_maf <- data.table::fread("GDCdata/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.csv")
#reads in as maf
maf <- read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

TP53mut <- subsetMaf(maf, gene = "TP53")
PIK3CAmut <- subsetMaf(maf, gene = "PIK3CA")
TTNmut <- subsetMaf(maf, gene = "TTN")

#maf@clinical.data$Tumor_Sample_Barcode
TP53mutbarcodes <- TP53mut@clinical.data$Tumor_Sample_Barcode #patients with mutated TP53
PIK3CAmutbarcodes <- PIK3CAmut@clinical.data$Tumor_Sample_Barcode
TTNmutbarcodes <- TTNmut@clinical.data$Tumor_Sample_Barcode

#barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
                     #"TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
                     #"TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
#load in HTSeq data
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts")
#GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

#Get counts for specific genes
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient #shorter barcodes to match those in MAF

#Access TP53 counts
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TP53_mask]

#TP53 counts for all counts in LumA barcodes and in TP53-mutated barcodes
LumA_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% LumAbarcodes %in% TP53mutbarcodes]
#TP53 counts for all counts in LumA barcodes but not in TP53-mutated barcodes
LumA_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% LumAbarcodes %in% !(colnames(htseq_counts) %in% TP53mutbarcodes)]

LumB_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% LumBbarcodes %in% TP53mutbarcodes]
LumB_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% LumBbarcodes %in% !(colnames(htseq_counts) %in% TP53mutbarcodes)]

Basal_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% Basalbarcodes %in% TP53mutbarcodes]
Basal_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% Basalbarcodes %in% !(colnames(htseq_counts) %in% TP53mutbarcodes)]

Her2_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% Her2barcodes %in% TP53mutbarcodes]
Her2_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% Her2barcodes %in% !(colnames(htseq_counts) %in% TP53mutbarcodes)]

Normal_mut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% Normalbarcodes %in% TP53mutbarcodes]
Normal_nonmut_TP53_counts <- htseq_counts[TP53_ENSG_name, colnames(htseq_counts) %in% Normalbarcodes %in% !(colnames(htseq_counts) %in% TP53mutbarcodes)]

#boxplot(mut_TP53_counts, nonmut_TP53_counts, names = c("mutated TP53 samples", "non-mutated TP53 samples"), 
        #ylab = "HTSeq mRNA counts for TP53 gene")
pdf("TP53_RNA_vs_mut_boxplot.pdf")
boxplot(LumA_mut_TP53_counts, LumA_nonmut_TP53_counts, LumB_mut_TP53_counts, LumB_nonmut_TP53_counts, Basal_mut_TP53_counts, 
        Basal_nonmut_TP53_counts, Her2_mut_TP53_counts, Her2_nonmut_TP53_counts, Normal_mut_TP53_counts, Normal_nonmut_TP53_counts,
        names = c("LumA mutated", "LumA non-mutated", "LumB mutated", "LumB non-mutated", "Basal mutated", "Basal non-mutated",
                  "Her2 mutated", "Her2 non-mutated", "Normal mutated", "Normal non-mutated"), 
        xlab = "Tumor subtype and TP53 mutation status", ylab = "HTSeq mRNA counts for TP53 gene")
dev.off()

#Access PIK3CA counts
#PIK3CA_mask <- rowData(sum_exp)$external_gene_name == "PIK3CA"
#PIK3CA_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[PIK3CA_mask]

#mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, colnames(htseq_counts) %in% PIK3CAmutbarcodes]
#nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, !(colnames(htseq_counts) %in% PIK3CAmutbarcodes)]

#mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, colnames(htseq_counts) %in% PIK3CAmutbarcodes]
#nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, !(colnames(htseq_counts) %in% PIK3CAmutbarcodes)]

#mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, colnames(htseq_counts) %in% PIK3CAmutbarcodes]
#nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, !(colnames(htseq_counts) %in% PIK3CAmutbarcodes)]

#mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, colnames(htseq_counts) %in% PIK3CAmutbarcodes]
#nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, !(colnames(htseq_counts) %in% PIK3CAmutbarcodes)]

#mut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, colnames(htseq_counts) %in% PIK3CAmutbarcodes]
#nonmut_PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name, !(colnames(htseq_counts) %in% PIK3CAmutbarcodes)]

#boxplot(mut_PIK3CA_counts, nonmut_PIK3CA_counts, names = c("mutated PIK3CA genes", "non-mutated PIK3CA genes"), 
        #ylab = "HTSeq mRNA counts for PIK3CA gene")



#Access TTN counts
#TTN_mask <- rowData(sum_exp)$external_gene_name == "TTN"
#TTN_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TTN_mask]

#mut_TTN_counts <- htseq_counts[TTN_ENSG_name, colnames(htseq_counts) %in% TTNmutbarcodes]
#nonmut_TTN_counts <- htseq_counts[TTN_ENSG_name, !(colnames(htseq_counts) %in% TTNmutbarcodes)]

#mut_TTN_counts <- htseq_counts[TTN_ENSG_name, colnames(htseq_counts) %in% TTNmutbarcodes]
#nonmut_TTN_counts <- htseq_counts[TTN_ENSG_name, !(colnames(htseq_counts) %in% TTNmutbarcodes)]

#mut_TTN_counts <- htseq_counts[TTN_ENSG_name, colnames(htseq_counts) %in% TTNmutbarcodes]
#nonmut_TTN_counts <- htseq_counts[TTN_ENSG_name, !(colnames(htseq_counts) %in% TTNmutbarcodes)]

#mut_TTN_counts <- htseq_counts[TTN_ENSG_name, colnames(htseq_counts) %in% TTNmutbarcodes]
#nonmut_TTN_counts <- htseq_counts[TTN_ENSG_name, !(colnames(htseq_counts) %in% TTNmutbarcodes)]

#mut_TTN_counts <- htseq_counts[TTN_ENSG_name, colnames(htseq_counts) %in% TTNmutbarcodes]
#nonmut_TTN_counts <- htseq_counts[TTN_ENSG_name, !(colnames(htseq_counts) %in% TTNmutbarcodes)]

#boxplot(mut_TTN_counts, nonmut_TTN_counts, names = c("mutated TTN genes", "non-mutated TTN genes"), 
        #ylab = "HTSeq mRNA counts for TTN gene")



