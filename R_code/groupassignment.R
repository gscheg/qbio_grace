# Question: What is the difference in survivability between old, mid, and young patients 
# with low vs. high gene expression of TP53, PIK3CA, MUC16?

# Installation
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if(!requireNamespace("devtools")) BiocManager::install(c("devtools"))
#if(!requireNamespace("robustbase"))BiocManager::install(c("robustbase"))
#library(devtools)
#library(robustbase)
#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
#if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
#if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
#if(!requireNamespace("arsenal"))install.packages(c("arsenal"))
#if(!requireNamespace("survival"))install.packages(c("survival"))
#if(!requireNamespace("survminer"))install.packages(c("survminer"))

# Loading packages
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(arsenal)
library(survival)
library(survminer)

# Barcodes
barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
                     "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
                     "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

# Accessing RNAseq data 
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  barcode = barcodes_rnaseq)
GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", barcode=barcodes_clinic, file.type = "xml")
GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

# Getting counts for specific genes, replacing column names with shorter barcodes (to match clinical)
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient

# Matching clinical data sample order to RNAseq sample order
row_order <- match(colnames(htseq_counts), clinic$bcr_patient_barcode)
clinic_ordered  <- clinic[row_order, ]
#clinic_matched$
# Get rid of nonmatching samples in clinical and htseq
matching <- which(clinic_ordered$bcr_patient_barcode %in% colnames(htseq_counts))
clinic_matched <- clinic_ordered[matching,]

# Accessing counts data for TP53, categorize expression, and add to clinical data
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ TP53_mask ]
TP53_counts <- htseq_counts[TP53_ENSG_name,clinic_matched$bcr_patient_barcode]
TP53_quartiles <- quantile(TP53_counts)
TP53_expression_level <- ifelse(TP53_counts > TP53_quartiles[4], "High", ifelse(TP53_counts < TP53_quartiles[2], "Low", "Mid"))
clinic_matched$TP53_counts <- TP53_counts
clinic_matched$TP53_expression = TP53_expression_level

# Accessing counts data for PIK3CA, categorize expression, and add to clinical data
PIK3CA_mask <- rowData(sum_exp)$external_gene_name == "PIK3CA"
PIK3CA_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ PIK3CA_mask ]
PIK3CA_counts <- htseq_counts[PIK3CA_ENSG_name,clinic_matched$bcr_patient_barcode]
PIK3CA_quartiles <- quantile(PIK3CA_counts)
PIK3CA_expression_level <- ifelse(PIK3CA_counts > PIK3CA_quartiles[4], "High", ifelse(PIK3CA_counts < PIK3CA_quartiles[2], "Low", "Mid"))
clinic_matched$PIK3CA_counts <- PIK3CA_counts
clinic_matched$PIK3CA_expression = PIK3CA_expression_level

# Accessing counts data for MUC16, categorize expression, and add to clinical data
MUC16_mask <- rowData(sum_exp)$external_gene_name == "MUC16"
MUC16_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ MUC16_mask ]
MUC16_counts <- htseq_counts[MUC16_ENSG_name,clinic_matched$bcr_patient_barcode]
MUC16_quartiles <- quantile(MUC16_counts)
MUC16_expression_level <- ifelse(MUC16_counts > MUC16_quartiles[4], "High", ifelse(MUC16_counts < MUC16_quartiles[2], "Low", "Mid"))
clinic_matched$MUC16_counts <- MUC16_counts
clinic_matched$MUC16_expression = MUC16_expression_level

# Adding age to clinical data
age_clinical = clinic_matched$age_at_initial_pathologic_diagnosis
clinic_matched$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

#Boxplots for gene expression of each gene as a function of age
boxplot(clinic_matched$PIK3CA_counts ~ clinic_matched$age_category, xlab = "age categories", ylab = "PIK3CA counts", 
        main = "PIK3CA mRNA copy counts distributed by age")
boxplot(clinic_matched$MUC16_counts ~ clinic_matched$age_category, xlab = "age categories", ylab = "MUC16 counts", 
        main = "MUC16 mRNA copy counts distributed by age")
boxplot(clinic_matched$TP53_counts ~ clinic_matched$age_category, xlab = "age categories", ylab = "TP53 counts", 
        main = "TP53 mRNA copy counts distributed by age")

#Attempt to make KM survival plots to display different survivabilities for each age group at a particular level of expression for a gene,
#but not yet successful
clinic_PIK3CAlow <- clinic_matched[PIK3CA_expression_level == "Low",]
clinic_matched[PIK3CA_expression_level == "Mid",]
clinic_matched[PIK3CA_expression_level == "High",]

TCGAanalyze_survival(clinic_matched[PIK3CA_expression_level == "Low",], clinic_PIK3CAlow$age_clinical, filename="./PIK3CAkaplanmeierbyage.pdf")
