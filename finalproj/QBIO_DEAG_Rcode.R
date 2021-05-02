if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
if (!require(maftools)) BiocManager::install("maftools")
library(maftools)
library(SummarizedExperiment)

#------------Recreate the Top 10 mutated genes MAF Summary Plot (so we can pick the top 2 mutated genes to investigate)------------#
mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="varscan2") # only need to query once
maf_dataframe = read.maf(mutation)
#See here: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ for information on each column

#to save the pdf, execute all your plotting code after the pdf() line:
pdf("final_project_maf_summary.pdf")
plotmafSummary(maf = maf_dataframe, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
# top 2 mutated genes are PIK3CA and TP53, and TTN

#------------Daniel's code

#install.packages("survival")
#install.packages("survminer")
#install.packages("arsenal")

#if you get an error saying you need a CRAN mirror, it is with the above. Opening an R terminal
# and manually installing the packages worked for me

library(survival)
library(survminer)
library(arsenal)   


clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
#GDCdownload( clin_query ) #should only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")  #these download clinical data into clinic matrix
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"   # fixes missing
# underscore, would interupt TCGA_analyze looking for that column otherwise


age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical <40, "Young", ifelse(age_clinical>=60, "Old", "Mid"))
#adds column classifying each patient as young, mid or old.

subtypes <- TCGAquery_subtype(tumor = "BRCA")
#grabs the database with PAM50 subtypes, clinic doesn't actually have

barcodes_subtypes = subtypes$patient
PAM50 = rep("NA",length(clinic$bcr_patient_barcode))
barcodes_clinic = clinic$bcr_patient_barcode
counter2=0
for (barcode1 in barcodes_clinic){
  counter = 0
  counter2= counter2 +1
  for (barcode2 in barcodes_subtypes){
    counter= counter +1
    if( barcode1==barcode2){
      PAM50[counter2] = subtypes$BRCA_Subtype_PAM50[counter]
    }
    
  }
}

clinic$PAM50 = PAM50

#this grabs the subtype column and puts it in clinic. Doing TCGA_analyze from subtypes gives errors that 
#I didn't want to try deciphering
#the nested for loops are needed cause subtypes is smaller. Doesn't have all the same patients, so need 
#verify barcodes


NA_list = c()

young_list = c()
counter = 0
for (tumor in PAM50){
  counter= counter +1
  if (tumor =="NA"){
    NA_list =c(NA_list,counter)
  }
}
#NA list now holds row numbers of Patients with no listed subtype
counter = 0
age_groups = clinic$age_category
for (age in age_groups){
  counter= counter +1
  if (age =="Young"){
    young_list =c(young_list,counter)
  }
}
#young list now holds row numbers of patients that are young

copy_clinic = clinic[-NA_list,] #copy is clinic but NA subtype rows excluded
young_clinic = copy_clinic[young_list,] #young is copy_clinic but just young (also no NA) 
old_clinic = copy_clinic[-young_list,] #old is copy_clinic but excluding young (also no NA)

TCGAanalyze_survival( copy_clinic, "PAM50", filename="PAM50_survival.pdf")

TCGAanalyze_survival( young_clinic, "PAM50", filename="young_PAM50_survival.pdf")

TCGAanalyze_survival( old_clinic, "PAM50", filename="older_PAM50_survival.pdf")

#TCGA_analyze grabs several columns and provided column to make Kaplan Meier survival plots



jpeg("distributetest.jpg")
barplot(table(copy_clinic$PAM50), main = "Distribution of Subtypes for Patients in TCGA Database",
        xlab = "PAM50 Subtype", ylab ="Number of Tumors")
dev.off()

#barplot makes a bargraph. Here shows distribution of subtypes. Using modified clinics to focus on young or old


jpeg("youngdistributetest.jpg")
barplot(table(young_clinic$PAM50), main = "Subtype distribution for young patients(<40 years)",
        xlab = "PAM50 Subtype", ylab ="Number of Tumors")
dev.off()

jpeg("olderdistributetest.jpg")
barplot(table(old_clinic$PAM50), main = "Subtype distribution for older patients(>=40 years)",
        xlab = "PAM50 Subtype", ylab ="Number of Tumors")
dev.off()


#------------Grace's code
#boxplots of transcriptomics vs. mutation data

#load in TCGA-BRCA clinical data
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
#read in maf
maf <- read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

#subset of MAF data for each gene of interest
TP53mut <- subsetMaf(maf, gene = "TP53")
PIK3CAmut <- subsetMaf(maf, gene = "PIK3CA")

#barcodes of patients with mutated gene for both genes of interest
TP53mutbarcodes <- TP53mut@clinical.data$Tumor_Sample_Barcode #patients with mutated TP53
PIK3CAmutbarcodes <- PIK3CAmut@clinical.data$Tumor_Sample_Barcode

#load in HTSeq data
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts")
GDCdownload(query)
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
