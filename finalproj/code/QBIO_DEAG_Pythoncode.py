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

protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()


clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12

plt.figure(figsize=(8,6))
sns.countplot(data = clinical_data, x = 'PAM50', order = clinical_data['PAM50'].value_counts().index)
plt.title("Distribution of Subtypes in CPTAC-BRCA")
plt.ylabel("Number of Patients")
plt.xlabel("Subtype by PAM50")

clinical_data['Age_in_years'].value_counts(dropna=False)

plt.figure(figsize=(8,6))
sns.countplot(data = clinical_data.loc[clinical_data['Age_in_years'] < 40], x = 'PAM50', 
              order = clinical_data.loc[clinical_data['Age_in_years'] < 40]['PAM50'].value_counts().index,
             color = 'steelblue')
plt.title("Distribution of Subtypes for Young Patients (<40 years) in CPTAC-BRCA")
plt.ylabel("Number of Patients")
plt.xlabel("Subtype by PAM50")

plt.figure(figsize=(8,6))
sns.countplot(data = clinical_data.loc[clinical_data['Age_in_years'] >= 60], x = 'PAM50', 
              order = clinical_data.loc[clinical_data['Age_in_years'] >= 60]['PAM50'].value_counts().index,
             color = 'steelblue')
plt.title("Distribution of Subtypes for Old Patients (>=60 years) in CPTAC-BRCA")
plt.ylabel("Number of Patients")
plt.xlabel("Subtype by PAM50")



clinical_data.shape

young = clinical_data.loc[clinical_data['Age_in_years'] < 40]
old = clinical_data.loc[clinical_data['Age_in_years'] >= 60]

clinical_data.head()

protein_data.head()

protein_young = protein_data.loc[young.index]
protein_old = protein_data.loc[old.index]

protein_young.shape

protein_old.shape

protein_young.mean().sort_values(ascending=False)[0:5]

protein_old.mean().sort_values(ascending=False)[0:5]

protein_data_marked = protein_data.copy()
protein_data_marked['Group'] = np.nan
for i in young.index:
    protein_data_marked.at[i, 'Group'] = "Young"
for j in old.index:
    protein_data_marked.at[j, 'Group'] = "Old"

sns.boxplot(x=protein_data_marked["Group"], y=protein_data_marked["PIK3CA"])
plt.title("PIK3CA Expression in Young and Old Patients in CPTAC-BRCA")

sns.boxplot(x=protein_data_marked["Group"], y=protein_data_marked["TP53"])
plt.title("TP53 Expression in Young and Old Patients in CPTAC-BRCA")

protein_data_marked = protein_data_marked.sort_index()
clinical_data = clinical_data.sort_index()

protein_data_marked['PAM50'] = clinical_data['PAM50']


sns.boxplot(x=protein_data_marked["PAM50"], y=protein_data_marked["PIK3CA"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("PIK3CA Expression by Subtype (PAM50) in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")

sns.boxplot(x=protein_data_marked["PAM50"], y=protein_data_marked["TP53"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("TP53 Expression by Subtype (PAM50) in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")


#protein_data_marked.loc[protein_data_marked['Group'] == "Young"]

sns.boxplot(x=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["PAM50"], y=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["PIK3CA"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("PIK3CA Expression by Subtype (PAM50) for Young Patients in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")

sns.boxplot(x=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["PAM50"], y=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["TP53"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("TP53 Expression by Subtype (PAM50) for Young Patients in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")

clinical_data.shape


from sklearn.manifold import TSNE
import umap

protein_data = protein_data.sort_index()

from sklearn.impute import SimpleImputer
values = protein_data.values
imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
imputed_values = imputer.fit_transform(values)

protein_data_imp = pd.DataFrame(imputed_values, columns = protein_data.columns, index= protein_data.index)

tsne = TSNE(n_components=2, random_state=7, perplexity=35, n_iter=1000, learning_rate=200)
tsne_res = tsne.fit_transform(protein_data_imp)

tsne_df = pd.DataFrame(data = tsne_res, columns = ['Dimension 1', 'Dimension 2'])
tsne_df = tsne_df.set_index(protein_data_imp.index)

final_tsne_df = pd.concat([tsne_df, clinical_data[['PAM50', 'Age_in_years']]], axis = 1)
final_tsne_df = pd.concat([final_tsne_df, protein_data_marked[["Group"]]], axis = 1)

curr_hue = 'PAM50'
plt.figure(figsize=(6,6))
sns.scatterplot(x="Dimension 1", y="Dimension 2", hue = curr_hue, data=final_tsne_df)
plt.title("t-SNE Clustering of Protein Expression Data (Colored by Subtype)")


