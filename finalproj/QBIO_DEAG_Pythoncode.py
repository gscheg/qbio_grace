import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

#Alyssa's code
cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

#Remove the "multi" part of the dataframe (the dataframes are MultIndex pandas dataframes)
protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()

clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12


# clinical_data.to_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/clinical_data.csv")
# protein_data.to_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/protein_data.csv")
# rna_data.to_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/rna_data.csv")

# clinical_data_readin = pd.read_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/clinical_data.csv", index_col=1) #The index_col=0 creates the index (rownames) as the first column. 


#Check that the patients are in the same order in rna_data and protein_data.
#We are looking at the gene expression and protein information of EACH patient, so the data needs to be in pairs.
#rna_data.index are the patients (not rna_data.columns)
assert list(protein_data.index) == list(rna_data.index)


#Creates scatter plots presenting the Spearman correlation rho value
def spearman_plots (gene_name, gene_rna_data, gene_protein_data, output_file):
    '''make scatterplots showing correlation between transcriptomic and proteomic data'''

    rho, spear_pvalue = stats.spearmanr( gene_rna_data, gene_protein_data ) #using scipy

    plt.figure( figsize=(10,10) ) #using matplotlib.pyplot
    plt.scatter( gene_rna_data, gene_protein_data )

    #trend line
    #https://stackoverflow.com/questions/41635448/how-can-i-draw-scatter-trend-line-on-matplot-python-pandas
    z = np.polyfit(gene_rna_data, gene_protein_data, 1)
    p = np.poly1d(z)
    plt.plot(gene_rna_data,p(gene_rna_data),"r--")

    title = "rho: {} for {}".format(rho, gene_name)
    plt.title(title)

    plt.xlabel("RNA")
    plt.ylabel("Protein")

    pngfile = '{}.png'.format( output_file)

    plt.savefig( pngfile, bbox_inches="tight" ) #Use plt.savefig() when saving figure in script (use plt.show() in jupyter ntbk)


#Access the transcriptomic and proteomic information of the specific genes (PIK3CA and TP53)
rna_pik3ca = rna_data.loc[:, "PIK3CA"] #want all rows of the [specifc gene] column. 
protein_pik3ca = protein_data.loc[:, "PIK3CA"]
spearman_plots("PIK3CA", rna_pik3ca, protein_pik3ca, "pik3ca_spear")

rna_tp53 = rna_data.loc[:, "TP53"]
protein_tp53 = protein_data.loc[:, "TP53"]
spearman_plots("TP53", rna_tp53, protein_tp53, "tp53_spear")



#Eric's code
cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

#Remove the "multi" part of the dataframe (the dataframes are MultIndex pandas dataframes)
protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()

clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12



#Check that the patients are in the same order in rna_data and protein_data.
#We are looking at the gene expression and protein information of EACH patient, so the data needs to be in pairs.
#rna_data.index are the patients (not rna_data.columns)
assert list(protein_data.index) == list(rna_data.index)


#Creates scatter plots presenting the Spearman correlation rho value
def spearman_plots (gene_name, gene_rna_data, gene_protein_data, output_file):
    '''make scatterplots showing correlation between transcriptomic and proteomic data'''

    rho, spear_pvalue = stats.spearmanr( gene_rna_data, gene_protein_data ) #using scipy

    plt.figure( figsize=(10,10) ) #using matplotlib.pyplot
    plt.scatter( gene_rna_data, gene_protein_data )

    #trend line
    #https://stackoverflow.com/questions/41635448/how-can-i-draw-scatter-trend-line-on-matplot-python-pandas
    z = np.polyfit(gene_rna_data, gene_protein_data, 1)
    p = np.poly1d(z)
    plt.plot(gene_rna_data,p(gene_rna_data),"r--")

    title = "rho: {} for {}".format(rho, gene_name)
    plt.title(title)

    plt.xlabel("RNA")
    plt.ylabel("Protein")

    pngfile = '{}.png'.format( output_file)

    plt.savefig( pngfile, bbox_inches="tight" ) #Use plt.savefig() when saving figure in script (use plt.show() in jupyter ntbk)


#Access the transcriptomic and proteomic information of the specific genes (PIK3CA and TP53)
rna_pik3ca = rna_data.loc[:, "PIK3CA"] #want all rows of the [specifc gene] column. 
protein_pik3ca = protein_data.loc[:, "PIK3CA"]
spearman_plots("PIK3CA", rna_pik3ca, protein_pik3ca, "pik3ca_spear")

rna_tp53 = rna_data.loc[:, "TP53"]
protein_tp53 = protein_data.loc[:, "TP53"]
spearman_plots("TP53", rna_tp53, protein_tp53, "tp53_spear")

