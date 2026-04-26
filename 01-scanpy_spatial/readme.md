#  Module 1: Spatial Transcriptomics Analysis (Scanpy)

---

##  Objective

This module follows the official Scanpy spatial transcriptomics workflow using 10x Genomics Visium data.

We aim to understand how gene expression varies across tissue space and how it can be used to identify biological structure.

### 📌 What you will learn:
- How to load spatial transcriptomics data (Visium)
- How to perform quality control (QC)
- How to normalize and preprocess gene expression data
- How to reduce dimensionality (PCA, UMAP)
- How to cluster tissue spots (Leiden algorithm)
- How to map gene expression back onto tissue images

---

##  Dataset

We use a 10x Genomics Visium dataset from a human lymph node.

This dataset contains:
- RNA expression matrix (genes × spots)
- Spatial coordinates for each spot
- Histology image (H&E stained tissue section)
- Metadata and QC statistics

###  Why this matters:
Visium data connects **gene expression + physical tissue location**, enabling spatial biology analysis.

---

##  Workflow Overview

```text
Data Loading → QC → Filtering → Normalization → HVGs → PCA → Neighbors → UMAP → Leiden → Visualization → Marker Genes
```

###  Explanation:
This pipeline converts raw gene counts into biologically meaningful clusters and spatial maps.

---

##  QC & Filtering

Quality control ensures that only high-quality tissue spots are used.

<img width="635" height="173" alt="gene by count scanpy" src="https://github.com/user-attachments/assets/d2152033-6aae-4158-995e-a54ab8dcdd50" />


###  What we check:
- Total RNA counts per spot
- Number of genes detected per spot
- Percentage of mitochondrial genes (cell stress indicator)

```python
import seaborn as sns

# QC distributions
sns.histplot(adata.obs["total_counts"])
sns.histplot(adata.obs["n_genes_by_counts"])

# Filtering low-quality spots
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)

# Remove stressed cells (high mitochondrial content)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

# Remove rarely expressed genes
sc.pp.filter_genes(adata, min_cells=10)
```

###  Why this matters:
Low-quality spots can distort clustering and lead to incorrect biological interpretation.

---

## ⚙️ Normalization & HVGs

Before analysis, data must be standardized.

###  Steps:
- Normalize counts per cell
- Log-transform data
- Select highly variable genes (HVGs)

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000
)
```

###  Why this matters:
- Removes sequencing depth bias
- Focuses on biologically informative genes
- Reduces noise for downstream analysis

---

##  PCA → Neighbors → UMAP → Leiden

This is the core analytical pipeline.

```python
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(
    adata,
    key_added="clusters",
    flavor="igraph",
    directed=False
)
```

###  Explanation of each step:

- **PCA** → compress gene expression into main variation axes
- **Neighbors graph** → connects similar spots
- **UMAP** → visualizes structure in 2D
- **Leiden clustering** → groups biologically similar regions

---

##  UMAP Visualizations

###  What UMAP shows:
A 2D map of transcriptionally similar tissue spots.

```python
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts"])
sc.pl.umap(adata, color="clusters")
```

<img width="709" height="185" alt="umap,scanpy" src="https://github.com/user-attachments/assets/922fa179-efd1-4dab-858c-baab848a2baf" />




###  Interpretation:
- Each point = a tissue spot
- Distance = similarity in gene expression
- Colors = discovered biological clusters

---

##  Spatial Visualization

Now we map results back onto real tissue.

<img width="640" height="296" alt="1st_spatial_scanpy" src="https://github.com/user-attachments/assets/8e1fa56f-ec5c-47c1-a17a-88b495a2e9d1" />

<img width="203" height="290" alt="2nd _spatial_scanpy" src="https://github.com/user-attachments/assets/bc98e0e0-d5d2-4b60-8e21-8314bfeab9c4" /><img width="312" height="293" alt="spatial_cluster scanpy" src="https://github.com/user-attachments/assets/1d67481b-8fa1-4367-975f-f2fe48250020" />


###  Explanation:
- Overlays clusters on histology image
- Shows where gene expression patterns exist in tissue
- Connects molecular data with physical structure

---

## Cluster marker genes

We visualize individual genes across tissue.

```python
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="9", n_genes=10, groupby="clusters")
)
```

<img width="296" height="369" alt="cluster marker genes_scanpy" src="https://github.com/user-attachments/assets/a3976ae2-0e00-4ab2-9a2e-31c7333bb4e0" />

### 📌 Interpretation:
- Shows where specific genes are active
- Helps identify functional tissue regions
- Confirms biological relevance of clusters

<img width="619" height="272" alt="marker_spatial cluter1 scanpy" src="https://github.com/user-attachments/assets/55dd1024-573a-4986-a0b8-4e7b08f1c228" />

<img width="597" height="267" alt="marker_spatial2" src="https://github.com/user-attachments/assets/c50f3a32-69a2-4e89-9ba7-37fd317313b0" />

---

## 🧠 Key Insights

- Spatial transcriptomics combines **gene expression + spatial location**
- Clusters represent real biological structures in tissue
- UMAP shows transcriptional similarity, not physical layout
- Spatial plots validate clustering results in real tissue context

---

## 🎯 Final Takeaway

Spatial transcriptomics allows us to understand:

> How gene activity is organized across physical tissue space, revealing hidden biological structure that cannot be seen with traditional RNA-seq.
