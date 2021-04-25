import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# cptac.download(dataset="Brca")
# br = cptac.Brca()

# protein_data = br.get_proteomics()

# # The dataframes are MultIndex pandas dataframes.
# # However, to teach the basics of pandas, we will remove the "multi" part of the dataframe.
# protein_data = protein_data.droplevel(1, axis=1)

# rna_data = br.get_transcriptomics()
# clinical_data = br.get_clinical()
# clinical_data["Age_in_years"] = clinical_data["Age.in.Month"] // 12  # fill in with correct math

# # A note about saving pandas

# clinical_data.to_csv("data/cptac_clinical.csv")
# protein_data.to_csv("data/cptac_protein.csv")
# rna_data.to_csv("data/cptac_rna.csv")

# Repeat with protein_data and rna_data if desired
# We would recommend creating a "data" folder where you can keep all these files organized

# Read back in your data
# Try with different arguments. The index_col=0 creates the index (rownames) as the first column.
# Double check that the index and columns are what you want
clinical_data = pd.read_csv("data/cptac_clinical.csv", index_col=0)
protein_data = pd.read_csv("data/cptac_protein.csv", index_col=0)
rna_data = pd.read_csv("data/cptac_rna.csv", index_col=0)

# Spearman Tutorial
# First we are going to look at individual genes
# Step 1: Check that the patients are in the same order.
# This is important because we are looking at the gene expression and protein information of EACH patient
# The data needs to be in pairs. This is easiest when the patients are in the same order in rna_data and protein_data

# Are the index or columns (rna_data.index or rna_data.columns) the patients?

# Write an assert statement that checks the patients rna_data are equal to the patients of protein_data
#Fill in parantheses
assert list(rna_data.index) == list(protein_data.index)


# Step 2: Access the transcriptomic and proteomic information of the specific gene (ESR1)

# Fill in the brackets []. We want ALL rows of the ESR1 column.
# Remember [row,col] format and : refers to ALL.
rna_esr1 = rna_data.loc[:, "ESR1"]
protein_esr1 = protein_data.loc[:, "ESR1"]

rho, spear_pvalue = stats.spearmanr(rna_esr1, protein_esr1)

# The order of transcriptomics and proteomics doesn't matter. Why?
rho_check, spear_pvalue_check = stats.spearmanr(protein_esr1, rna_esr1)

assert rho == rho_check

# Step 3: Plot the data
# You can present the rho value in various ways. Table, scatterplot, etc.
# Here you will make a scatterplot

plt.figure(figsize=(10, 10))

# Replace x and y with appropriate variables
plt.scatter(rna_esr1, protein_esr1)

title = "rho: {} for ESR1".format(rho)  # This is string formatting. The variable in the () will print in the {}
plt.title(title)

# Fill in informative x and y labels
plt.xlabel("RNA expression levels")
plt.ylabel("Protein expression levels")

# plt.show()  # Comment out when running in script
plt.savefig("figs/scatter_cors.png", bbox_inches="tight")  # Use this when saving figure in script


# Repeat with only young and only old patients for practice

# Access only young and only old patients via a boolean mask

# What column of clinical_data is referring to age?
young_mask = clinical_data["Age_in_years"] < 40.0
old_mask = clinical_data["Age_in_years"] > 60.0

# Check for understanding: Why do the below lines work?
rna_esr1_young = rna_data["ESR1"][young_mask]
protein_esr1_young = protein_data["ESR1"][young_mask]

# We want all patients of the ESR1 column
rna_esr1_old = rna_data["ESR1"][old_mask]  # fill in here
protein_esr1_old = protein_data["ESR1"][old_mask]  # fill in here

fig, ax = plt.subplots(2, constrained_layout=True)

ax[0].scatter(rna_esr1_young, protein_esr1_young)
ax[1].scatter(rna_esr1_old, protein_esr1_old)

rho_young, spear_pvalue_young = stats.spearmanr(rna_esr1_young, protein_esr1_young)
rho_old, spear_pvalue_old = stats.spearmanr(rna_esr1_old, protein_esr1_old)
rho_check_young, spear_pvalue_check_young = stats.spearmanr(protein_esr1_young, rna_esr1_young)
rho_check_old, spear_pvalue_check_old = stats.spearmanr(protein_esr1_old, rna_esr1_old)

assert rho_young == rho_check_young
assert rho_old == rho_check_old

ax[0].set(title=f"Young, rho = {rho_young}", xlabel="RNA", ylabel="Protein")
ax[1].set(title=f"Old, rho = {rho_old}", xlabel="RNA", ylabel="Protein")

# plt.show()
plt.savefig("figs/scatter_cors_yo.png", bbox_inches="tight")

# Create the two scatter plots!

# Spearman Heatmaps
# Heatmap: Transcriptomics vs. Proteomics

# Step 2: Find the genes that are in the both protein_data and rna_data

# There are NA's in the protein_data and rna_data
# We will just drop any genes (columns) that have NA's
protein_data.dropna(inplace=True, axis=1)
rna_data.dropna(inplace=True, axis=1)

# Notice that the there are repeated genes in the protein data
# Does this make sense?
print(protein_data.shape[1])  # number of columns (genes)
print(len(np.unique(protein_data.columns)))  # number of UNIQUE columns (genes)

unique_genes = np.intersect1d(rna_data.columns, protein_data.columns)

# Fill in the []. Use : and unique_genes.
protein_data_genes = protein_data.loc[:, unique_genes]

# Check your undertanding of the following lines!
genes_with_repeats = list(protein_data_genes.columns)
rna_data_genes = rna_data.loc[:, genes_with_repeats]

# Check for understanding: What do the assert statements check?
assert list(rna_data_genes.index) == list(protein_data_genes.index)
assert list(rna_data_genes.columns) == list(protein_data_genes.columns)

# Step 3: Calculate the Spearman correlations

# We need to calculate a spearman rho for each pair of genes (transcriptomics and proteomics)
# this calls for nested for loops!

# What are the dimensions of your final graph?
genes_with_repeats_subset = genes_with_repeats[1:21]
n = len(genes_with_repeats_subset)  # fill in with the dimension of one axis. The graph is square.
final_graph = np.zeros((n, n))

for i, geneA in enumerate(genes_with_repeats_subset):
    # we will fill in row by row
    graph_row = []
    for j, geneB in enumerate(genes_with_repeats_subset):

        rho, spear_pvalue = stats.spearmanr(rna_data_genes.iloc[:, i], protein_data_genes.iloc[:, j])

        # append rho to graph_row here
        graph_row.append(rho)

    # set row i of final_graph equal to graph_row
    # remember [row, col] format, with : referring to all rows or all cols
    final_graph[i, :] = graph_row

# Step 4: Create the Heatmap

plt.figure(figsize=(10, 10))
ax_spear = sns.heatmap(final_graph, cmap='RdBu', xticklabels=genes_with_repeats_subset,
                       yticklabels=genes_with_repeats_subset, center=0, vmin=-1, vmax=1)
ax_spear.set_title("Spearman Rho Values")
ax_spear.set_ylabel("Genes Y")
ax_spear.set_xlabel("Genes X")

# plt.show()  # Use in Jupyter, Comment out if running full script
plt.savefig("figs/heatmap.png", bbox_inches="tight") #Use to save figure if running full script

ax = sns.clustermap(final_graph)

reorder_rows, reorder_cols = ax.dendrogram_row.reordered_ind, ax.dendrogram_col.reordered_ind

annot_x = [genes_with_repeats_subset[x] for x in reorder_cols]
annot_y = [genes_with_repeats_subset[y] for y in reorder_rows]


# Distribution of Spearman
# We already set up the rna_data and protein_data
assert list(rna_data_genes.index) == list(protein_data_genes.index)
assert list(rna_data_genes.columns) == list(protein_data_genes.columns)

# We only need ONE for loop this time... Why?

rho_list = []
for ii, gene in enumerate(list(rna_data_genes)):
    rho, spear_pvalue = stats.spearmanr(rna_data_genes.iloc[:, ii], protein_data_genes.iloc[:, ii])
    # append rho to rho_list
    rho_list.append(rho)


plt.figure(figsize=(10, 10))
plt.hist(rho_list, bins=20)
plt.xlim(-1, 1)
plt.xlabel("rho value")
plt.ylabel("Number of Genes with rho")

# plt.show()  # Use in Jupyter, Comment out if running full script
plt.savefig("figs/hist.png", bbox_inches="tight")  # Use to save figure if running full script
