#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#data download from GEO 
# Euctopic and ectopic human endometrium (endometriosis)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11691


# In[31]:


import GEOparse

# Replace 'GSE11691' with your GEO series ID
gse_id = 'GSE11691'
gse = GEOparse.get_GEO(geo=gse_id, destdir="./dbGSEdata")


# In[32]:


import GEOparse

# Load the GEO dataset
gse = GEOparse.get_GEO(filepath="./dbGSEdata/GSE11691_family.soft.gz")


# In[39]:


import pandas as pd

# Extract samples and metadata
samples = gse.gsms
metadata = pd.DataFrame({
    sample: gse.gsms[sample].metadata
    for sample in gse.gsms.keys()
}).T

# Extract expression data (assuming it is stored in the table_rows attribute)
expression_data = pd.DataFrame({
    gsm: gse.gsms[gsm].table.set_index('ID_REF')['VALUE']
    for gsm in gse.gsms.keys()
})

print(expression_data)


# In[36]:


# Check for missing values
missing_genes = expression_data.index[expression_data.isnull().any(axis=1)]
print("Missing genes:", missing_genes)


# In[37]:


# Option 1: Remove missing values
expression_data_clean = expression_data.dropna()
# print(expression_data_clean)
# Option 2: Fill missing values (e.g., with zero)
# expression_data_clean = expression_data.fillna(0)


# In[40]:



import scanpy as sc

# Convert to AnnData
adata = sc.AnnData(X=expression_data_clean.T)

# Proceed with analysis
# Quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Filter the data
adata = adata[:, adata.var.highly_variable]

# Scale the data
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Compute UMAP
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=10)
sc.tl.umap(adata)


# In[41]:


# Use available gene that is present
available_genes = adata.var_names[:5]  # List the first five available genes
print("Available genes for plotting:", available_genes)

# Plot using an available gene
sc.pl.umap(adata, color=available_genes[0])  # Use the first available gene


# In[42]:


# Plot clusters instead of specific genes if needed

# Cluster the data
sc.tl.leiden(adata, resolution=0.5)

# Plot clusters
sc.pl.umap(adata, color='leiden')

# Find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# Plot the top markers
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)


# In[56]:


# K-means Clustering
from sklearn.cluster import KMeans

# Run K-means clustering
kmeans = KMeans(n_clusters=4, random_state=0).fit(adata.X)

# Add K-means clusters to adata
adata.obs['kmeans'] = kmeans.labels_.astype(str)

# Plot clusters
sc.pl.umap(adata, color='kmeans')

# Find marker genes for K-means clusters
sc.tl.rank_genes_groups(adata, 'kmeans', method='t-test')

# Plot the top markers
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)


# In[52]:


# adjusting the number of clusters or removing clusters with too few samples before ranking genes.

# Run K-means clustering with fewer clusters
kmeans = KMeans(n_clusters=4, random_state=0).fit(adata.X)

# Add K-means clusters to adata
adata.obs['kmeans'] = kmeans.labels_.astype(str)

# Plot clusters
sc.pl.umap(adata, color='kmeans')

# Check cluster sizes
print(adata.obs['kmeans'].value_counts())


# In[53]:


# Filter out clusters with fewer than a threshold (e.g., 5 samples)
valid_clusters = adata.obs['kmeans'].value_counts()[adata.obs['kmeans'].value_counts() > 1].index
adata_filtered = adata[adata.obs['kmeans'].isin(valid_clusters)].copy()

# Find marker genes for valid clusters
sc.tl.rank_genes_groups(adata_filtered, 'kmeans', method='t-test')

# Plot the top markers
sc.pl.rank_genes_groups(adata_filtered, n_genes=20, sharey=False)


# In[64]:


# Hierarchical Clustering
# Define a threshold for minimum cluster size
min_cluster_size = 4

# Filter clusters with enough samples
valid_clusters = adata.obs['hierarchical'].value_counts()[
    adata.obs['hierarchical'].value_counts() >= min_cluster_size
].index

# Subset the data
adata_filtered = adata[adata.obs['hierarchical'].isin(valid_clusters)].copy()

# Convert to categorical
adata_filtered.obs['hierarchical'] = adata_filtered.obs['hierarchical'].astype('category')

# Plot clusters
sc.pl.umap(adata_filtered, color='hierarchical')

# Find marker genes for valid clusters
sc.tl.rank_genes_groups(adata_filtered, 'hierarchical', method='t-test')

# Plot the top markers
sc.pl.rank_genes_groups(adata_filtered, n_genes=20, sharey=False)



# In[ ]:




