# Installation of the required python libraries
!pip install scanpy
!pip install anndata
!pip3 install igraph
!pip install celltypist
!pip install decoupler
!pip install fa2-modified

# Import core single cell tools

import scanpy as sc
import anndata as ad

# Dataset (which was uploaded on Google Colab session)
!wget /content/bone_marrow.h5ad
bone_marrow_adata = sc.read_h5ad("/content/bone_marrow.h5ad")
print(bone_marrow_adata)


# the dimensions of the dataset
bone_marrow_adata.shape
bone_marrow_adata.obs

# look at the first 5 rows describing the genes in the dataset
bone_marrow_adata.var.head()

# look at the first 5 rows describing the cells (ID) in the dataset
bone_marrow_adata.obs.head()

# proper dataframe format
bone_marrow_adata.to_df()

#14783 rows × 17374 columns
# Meaning 14783 cells 17374 genes


#QUALITY CONTROL (QC)

# Unique names to avoid duplicate entries
bone_marrow_adata.var_names_make_unique()
bone_marrow_adata.obs_names_make_unique()

# Search for possible contamination from dying cells (Mitochondrial (MT) content), ribosomal transcripts(RIBO) or hemoglobin(HB)
# The thresholds are : (MT < 5%, RB <10%, and HB < 5%)
bone_marrow_adata.var['MT'] = bone_marrow_adata.var_names.str.startswith("MT-")
bone_marrow_adata.var['RIBO'] = bone_marrow_adata.var_names.str.startswith("RPS", "RPL")
bone_marrow_adata.var['HB'] = bone_marrow_adata.var_names.str.startswith("^HB[^(P)]")

# Calculate the qc metrics
sc.pp.calculate_qc_metrics(
    bone_marrow_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)

# Check if qc metrics is included in the headers of obs (cell dataset)
bone_marrow_adata.obs.head()

# Check if qc metric is included in the headers of var (gene dataset)
bone_marrow_adata.var.head()


# Check results for dying cells(MT content)
mt_genes = bone_marrow_adata.var[bone_marrow_adata.var['MT']]
mt_genes

print("First 5 rows of 'pct_counts_MT':")
print(bone_marrow_adata.obs['pct_counts_MT'].head())

print("\nDescriptive statistics for 'pct_counts_MT':")
print(bone_marrow_adata.obs['pct_counts_MT'].describe())

# Check results for hemoglobin contamination (HB content). No HB contamination
hb_genes = bone_marrow_adata.var[bone_marrow_adata.var['HB']]
hb_genes

print("First 5 rows of 'pct_counts_HB':")
print(bone_marrow_adata.obs['pct_counts_HB'].head())

print("\nDescriptive statistics for 'pct_counts_HB':")
print(bone_marrow_adata.obs['pct_counts_HB'].describe())

# Check results for ribosomal transcripts.
ribo_genes = bone_marrow_adata.var[bone_marrow_adata.var['RIBO']]
ribo_genes  

print("First 5 rows of 'pct_counts_RIBO':")
print(bone_marrow_adata.obs['pct_counts_RIBO'].head())

print("\nDescriptive statistics for 'pct_counts_RIBO':")
print(bone_marrow_adata.obs['pct_counts_RIBO'].describe())


# Visualise the each aspect of the data through violin plots

# Violin plot to check the no. of genes expressed in each cell
sc.pl.violin(
    bone_marrow_adata,
    ["n_genes_by_counts"],
    jitter=0.4,
    multi_panel=False,
)

# Violin plot to see the total no. of Unique Molecular Identifiers (UMIs) detected in a cell
sc.pl.violin(
    bone_marrow_adata,
    ["total_counts"],
    jitter=0.4,
    multi_panel=False,
)	

# Violin plot to check the level of mitochondrial (MT)content
sc.pl.violin(
    bone_marrow_adata,
    ["pct_counts_MT"],
    jitter=0.4,
    multi_panel=False,
)

# Violin plot to check the level of ribosomal(RIBO) transcripts 
sc.pl.violin(
    bone_marrow_adata,
    ["pct_counts_RIBO"],
    jitter=0.4,
    multi_panel=False,
)

# Violin plot to check the level of hemoglobin contamination (HB)
sc.pl.violin(
    bone_marrow_adata,
    ["pct_counts_HB"],
    jitter=0.4,
    multi_panel=False,
)

# No filtration step is required as the levels of MT, RIBO and HB contaminations is below the desired thresholds for each (MT < 5%, RB <10%, and HB < 5%)

# Doublet detection and removal
sc.pp.scrublet(bone_marrow_adata)

# NORMALISATION

# Save a copy of the data
bone_marrow_adata.layers["counts"] = bone_marrow_adata.X.copy()

# Normalise the data
sc.pp.normalize_total(bone_marrow_adata)

# Logarithmize the data
sc.pp.log1p(bone_marrow_adata)

# Feature selection: selecting the top 1000 most variable genes
sc.pp.highly_variable_genes(bone_marrow_adata, n_top_genes=1000)

# Visualization of the normalisation results: after (left) and before (right) normalization
# left is normalized
# right is not
sc.pl.highly_variable_genes(bone_marrow_adata ) 

# DIMENTIONALITY REDUCTION & CLUSTERING

# Performing principle component analysis of the dataset
sc.tl.pca(bone_marrow_adata)

# Visualise the identified PCAs (top 10 PCAs)
sc.pl.pca_variance_ratio(bone_marrow_adata, n_pcs=10, log=False)


# Use PCA cooardinates and calculate proximity of the data points between each other for further umap dimensionality reduction.
sc.pp.neighbors(bone_marrow_adata)
sc.tl.umap(bone_marrow_adata)

# Visualise the umap created in the previous step (Scatter plot)
sc.pl.umap(bone_marrow_adata)

# Graph-based clustering using leiden algorithm.
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2)

# Visualization of the leiden clustering.
sc.pl.umap(
    bone_marrow_adata,
    color=["leiden"],
    size=8,
)

# For better visualisation of the previous plot
sc.pl.umap(
    bone_marrow_adata,
    color=["leiden"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    ncols = 1
)

# Resolution adjustment for Leiden clustering
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_02", resolution=0.02)
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res2", resolution=2)

# Visualization of the clustering with all 3 resolutions

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res0_02"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1
)

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res0_5"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1
)

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res0_5"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1,
    legend_loc="on data"
)

# Based on the umaps with all three resolutions, choose the best resolution for further analysis. The resolution chosen in this analysis is 0.5 (res0_5)

# CELL ANNOTATION AND CLUSTER ANNOTATION

# steps needed to translate Ensemble IDs to gene names so that the Decoupler tool can identify them and moves the analysis further ahead

!wget wget -O result.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>'

import pandas as pd

ensembl_var = pd.read_csv('/content/result.txt', header = None)

ensembl_var.columns = ['ensembl_gene_id', 'gene_name']

ensembl_var.head(3)

import decoupler as dc

# Query Omnipath and get PanglaoDB
markers = dc.op.resource(name="PanglaoDB", organism="human")

# Keep canonical cell type markers alone
#markers = markers[markers["canonical_marker"]]

# Remove duplicated entries
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

# Format because dc only accepts cell_type and genesymbol

markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target"]]


markers.head()

# Correct target to ensemble
markers = markers.merge(ensembl_var, left_on="target", right_on="gene_name", how="left")
markers = markers.drop(columns=["target"])

# Remove duplicated entries
markers = markers[~markers.duplicated(["source", "ensembl_gene_id"])]

# Format because dc only accepts cell_type and genesymbol
markers = markers.rename(columns={"source": "source", "ensembl_gene_id": "target"})

markers = markers[["source", "target"]]
markers = markers.dropna()

markers.head()

# Load the gene expression matrix into dc
dc.mt.ulm(data=bone_marrow_adata,
          net=markers,
          tmin = 3)

# Retrieve the score for each cell type
score = dc.pp.get_obsm(bone_marrow_adata, key="score_ulm")
score

# Preview the data
bone_marrow_adata.obsm["score_ulm"].head()
bone_marrow_adata.obsm["score_ulm"].columns

# Create the UMAP plot with resolution equal to 0.5 (res0_5)
score.obs['leiden_res0_5'] = bone_marrow_adata.obs['leiden_res0_5']
sc.pl.umap(score, color=["Neutrophils", "leiden_res0_5"], cmap="RdBu_r")


# Create object score and and copy clustering labels in it.
# Copy the leiden_res0_5 column from bone_marrow_adata.obs to score.obs
score.obs['leiden_res0_5'] = bone_marrow_adata.obs['leiden_res0_5']

# Rank genes
# Find differential abundance of cell types between clusters.
bone_marrow_adata_rank = dc.tl.rankby_group(score, groupby="leiden_res0_5", reference="rest", method="t-test_overestim_var")

# Filter the results of the differentially abundant cell types to only include cell types that have a positive test statistic.
bone_marrow_adata_rank = bone_marrow_adata_rank[bone_marrow_adata_rank["stat"] > 0]

# Show top rows of the result.
bone_marrow_adata_rank.head()

# Create a mapping of the cluster annotation (labels).
cluster_annotations = bone_marrow_adata_rank[bone_marrow_adata_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()

cluster_annotations

# Create a new column called 'cell_type' in .obs
bone_marrow_adata.obs['cell_type'] = bone_marrow_adata.obs['leiden_res0_5'].map(cluster_annotations)

# UMAP 2D visualization of the clusters annotated by cell type.
sc.pl.umap(
    bone_marrow_adata,
    color=['cell_type'],
    legend_loc="on data"
)

# Create a subset for multiple genes in the 'source' column. Gene names of the clusters that have been found in the previous step.
available_genes = set(bone_marrow_adata.var_names)

#Neutrophils
neutro_markers = markers[markers['source'].isin(['Neutrophils'])]['target']
neutro_markers = neutro_markers[neutro_markers.isin(available_genes)]
neutro_markers = neutro_markers.drop_duplicates()

display(neutro_markers)

# Gamma delta T cells
gmtc_markers = markers[markers['source'].isin(['Gamma delta T cells'])]['target']
gmtc_markers = gmtc_markers[gmtc_markers.isin(available_genes)]
gmtc_markers = gmtc_markers.drop_duplicates()

display(gmtc_markers)

#T memory cells
tmem_markers = markers[markers['source'].isin(['T memory cells'])]['target']
tmem_markers = tmem_markers[tmem_markers.isin(available_genes)]
tmem_markers = tmem_markers.drop_duplicates()

display(tmem_markers)

# NK cells
nk_cell_markers = markers[markers['source'].isin(['NK cells'])]['target']
nk_cell_markers = nk_cell_markers[nk_cell_markers.isin(available_genes)]
nk_cell_markers = nk_cell_markers.drop_duplicates()

display(nk_cell_markers)

# Nuocytes
nuo_markers = markers[markers['source'].isin(['Nuocytes'])]['target']
nuo_markers = nuo_markers[nuo_markers.isin(available_genes)]
nuo_markers = nuo_markers.drop_duplicates()

display(nuo_markers)

# B cells naive
bcn_markers = markers[markers['source'].isin(['B cells naive'])]['target']
bcn_markers = bcn_markers[bcn_markers.isin(available_genes)]
bcn_markers = bcn_markers.drop_duplicates()


display(bcn_markers)

# Platelets
plt_markers = markers[markers['source'].isin(['Platelets'])]['target']
plt_markers = plt_markers[plt_markers.isin(available_genes)]
plt_markers = plt_markers.drop_duplicates()


display(plt_markers)

# Plasma cells
plasma_markers = markers[markers['source'].isin(['Plasma cells'])]['target']
plasma_markers = plasma_markers[plasma_markers.isin(available_genes)]
plasma_markers = plasma_markers.drop_duplicates()


display(plasma_markers)

# Monocytes
mono_markers = markers[markers['source'].isin(['Monocytes'])]['target']
mono_markers = mono_markers[mono_markers.isin(available_genes)]
mono_markers = mono_markers.drop_duplicates()

display(mono_markers)


# Another way to visualise the cell types
marker_genes_dict = {
    "Neutrophils": neutro_markers.head().tolist(),
    "Gamma delta T cells": gmtc_markers.head().tolist(),
    "T memory cells": tmem_markers.head().tolist(),
    "NK cells": nk_cell_markers.head().tolist(),
    "Nuocytes": nuo_markers.head().tolist(),
    "B cells naive": bcn_markers.head().tolist(),
    "Plasma cells": plasma_markers.head().tolist(),
    "Monocytes": mono_markers.head().tolist(),

}

# Visualise the marker_genes dictionary created in the previous step

# Dot plot
sc.pl.dotplot(bone_marrow_adata, marker_genes_dict, "cell_type", dendrogram=True)

# Stacked Violin Plot
sc.pl.stacked_violin(
    bone_marrow_adata, marker_genes_dict, groupby="leiden_res0_5",  dendrogram=True
)

# Matrix plot
sc.pl.matrixplot(
    bone_marrow_adata,
    marker_genes_dict,
    "leiden_res0_5",
    dendrogram=True,
    cmap="Blues",
)


# Heat map
sc.pl.heatmap(
    bone_marrow_adata, marker_genes_dict, groupby="leiden_res0_5", cmap="viridis", dendrogram=True
)



# Questions relating to Stage two task

1. What cell types did you identify?

'0': 'Neutrophils',
'1': 'Gamma delta T cells',
'2': 'T memory cells',
'3': 'NK cells',
'4': 'Nuocytes',
'5': 'B cells naive',
'6': 'Platelets',
'7': 'Plasma cells',
'8': 'Monocytes'

2. Explain the biological role of each cell type
 
 '0': Neutrophils: Most abundant cell type among leukocytes. Modulate inflammation, carries out phagocytosis, tissue repair and clearance. First responders to infection

 '1': Gamma delta T cells: First line of defense against microbes. Bridges innate and adaptive immunity. Anti-tumour activity. 

 '2': T memory cells: Provides long-term immunity. Rapid secondary immune response. Immune regulation and tissue surveillance.

 '3': NK cells: First line of defense against viruses. Tumour surveillance and tissue homeostasis. 

 '4': Nuocytes: Mediates type-2 immunity. Primarily involved in defense against helminth infections. Mediates allergic asthma and inflammation.

 '5': B cells naive: Carries out antigen recognition and surveillance. Antigen-presenting cells. Upon encounter with an antigen they activate and differentiate into plasma cells or memory cells

 '6': Platelets: Maintains homeostasis by forming clots (coagulation). Mediates wound repair and regeneration. 

 '7': Plasma cells: Main player of the humoral immune system. Antibody-producing cells. Carry out immunomodulation by producing anti-inflammatory cytokines (IL-10, IL-35)

 '8': Monocytes: Highly efficient phagocytes. Mediates regulation of inflammation. Carries out tumour surveillance and tissue homeostasis. 


3. Is the tissue source really bone marrow? Justify your answer
 
 No, the tissue source is not bone marrow.

Bone marrow is the site for haematopoiesis, i.e. the process of creating all mature blood cell types and also supporting cells called the stromal cells. Therefore bone marrow contains two primary populations of cell types : Haematopoietic stem cells (HSCs) which produce blood cells and Mesenchymal stem cells (MSCs) which produce the stromal cells. Since these two types of progenitor cells were absent from the provided single cell dataset, it can be said that the  tissue source was not bone marrow.


4. Based on the relative abundance of cell types, is the patient healthy or infected?

The UMAP visualisation of the annotated cell clusters displays a large number of cells clustering under the gamma delta T cells and T memory cells groups respectively. While T memory cells are involved in long term immunity, typically indicating a sign of past infection or vaccination, but a large number of T memory cells can also be a sign of ongoing infections in some cases. In healthy individuals gamma delta T cells normally range between 1 - 10%. Therefore, a large cell cluster under the gamma delta T cell group strongly suggests to the presence of an ongoing infection (mostly microbial). 
In addition to the above mentioned evidences, the subsequent dot plot created from the marker genes dictionary indicates a high expression of the gene ENSG00000277632 (which is the ensemble ID for the gene CCL3 - C-C motif chemokine ligand 3) by neutrophils. When a pathogen is detected, neutrophils rapidly produce and release large amounts of CCL3 to recruit other immune cells, such as monocytes and dendritic cells, to the site of infection. 

Therefore, based on the relative abundance of cell types and their gene expression patterns, it would be safe to say that the patient is infected. 















