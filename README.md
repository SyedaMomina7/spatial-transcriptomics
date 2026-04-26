# 🧬 Spatial Transcriptomics Complete Pipeline

**Scanpy · Squidpy · Visium · Xenium**

---

##  Overview

Spatial transcriptomics is a revolutionary technology that enables the measurement of gene expression while preserving the spatial location of cells within tissue.

Unlike traditional single-cell RNA sequencing (scRNA-seq), which loses positional information, spatial transcriptomics allows us to:

* Understand **tissue architecture**
* Identify **spatially distinct cell populations**
* Study **cell–cell interactions**
* Link **gene expression with histological structure**

This repository presents a **complete end-to-end spatial transcriptomics pipeline**, integrating multiple tools and datasets across four major workflows.

---

##  Objectives

* Perform spatial transcriptomics analysis using **Scanpy** and **Squidpy**
* Analyze **Visium datasets (fluorescence + H&E imaging)**
* Explore **high-resolution Xenium spatial data**
* Understand spatial clustering, neighborhood interactions, and tissue structure

---

##  Tools & Technologies

| Tool       | Purpose                                                  |
| ---------- | -------------------------------------------------------- |
| Scanpy     | Core preprocessing, clustering, dimensionality reduction |
| Squidpy    | Spatial analysis, graph construction, image integration  |
| AnnData    | Efficient data structure for omics data                  |
| Matplotlib | Visualization                                            |

---

## 🔬 Pipeline Overview

![Pipeline](assets/pipeline.png)

**Steps:**

1. Data Loading (Visium / Xenium)
2. Quality Control & Filtering
3. Normalization & Scaling
4. Dimensionality Reduction (PCA, UMAP)
5. Clustering
6. Spatial Visualization
7. Spatial Analysis (Neighborhood, Image Features)

---

## Key Workflows

### 1️⃣ Scanpy Spatial Basics

* Preprocessing and clustering
* UMAP visualization
* Spatial mapping of clusters

### 2️⃣ Visium Fluorescence (Squidpy)

* Spatial neighbor graph
* Fluorescence image integration
* Spatial statistics

### 3️⃣ Visium H&E Analysis

* Histology image integration
* Gene expression overlay on tissue
* Morphology vs expression

### 4️⃣ Xenium Analysis

* High-resolution spatial transcriptomics
* Cell-level / subcellular resolution
* Advanced spatial mapping

---

## 📊 Key Results & Visual Insights

### 🔹 1. UMAP Clustering (Global Gene Expression Structure)

👉 From **Scanpy module**

<img width="709" height="185" alt="umap,scanpy" src="https://github.com/user-attachments/assets/58fddc7f-a947-4971-af85-44f3ca10f3e1" />

📌 **Insight:**
Cells cluster based on gene expression profiles, revealing distinct biological populations.

---

### 🔹 2. Spatial Cluster Mapping on Tissue

👉 From **Visium H&E 
<img width="890" height="166" alt="spatial_cluster h e" src="https://github.com/user-attachments/assets/67d9ba2f-ca1a-4d6b-8953-01f53a1c4241" />


📌 **Insight:**
Clusters are not random — they align with **physical regions of the tissue**, showing spatial organization of biology.


### 🔹 3. High-Resolution Xenium Spatial Map

 From **Xenium module**



 **Insight:**
Unlike Visium (spot-based), Xenium provides **near single-cell resolution**, enabling precise spatial localization of gene expression.

---

##  Key Concepts Learned

* Spatial clustering of cells
* Spatially variable genes
* Neighborhood graph construction
* Integration of gene expression with tissue images
* High-resolution spatial transcriptomics

---

##  Repository Structure

```
## 📁 Repository Structure

spatial-transcriptomics/
├── 01-scanpy_spatial/
│   ├── visuals/
│   │   ├── 1st_spatial_scanpy.png
│   │   ├── 2nd_spatial_scanpy.png
│   │   ├── cluster_marker_genes_scanpy.png
│   │   ├── gene_by_count_scanpy.png
│   │   ├── marker_spatial_cluster1_scanpy.png
│   │   ├── marker_spatial2.png
│   │   ├── spatial_cluster_scanpy.png
│   │   └── umap_scanpy.png
│   ├── Analysis_and_visualization_of_spatial_transcriptomics_data.ipynb
│   └── readme.md
├── 02-Visium_fluorescence/
│   ├── visuals/
│   │   ├── Feature_Types_Comparison_flu.png
│   │   ├── Fluorescence_Channels.png
│   │   ├── features_vs_clusters_flu.png
│   │   ├── img_seg_flu.png
│   │   └── spatial_cluster_flu.png
│   ├── Analyze_Visium_fluorescence_data.ipynb
│   └── readme.md
├── 03-Visium_HnE/
│   ├── visuals/
│   │   ├── Co-occurrence_Heatmap_hne.png
│   │   ├── Gene_vs_Image_Clustering_hne.png
│   │   ├── moron.png
│   │   ├── spatial_cluster_hne.png
│   │   └── spatial_cluster_map_hne.png
│   ├── Analyze_Visium_HnE_data.ipynb
│   └── readme.md
└── 04-Xenium_data/
    ├── visuals/
    │   ├── QC_distributions_histogram_X.png
    │   ├── UMAP_Total_Counts_Leiden_Clusters_X.png
    │   ├── spatial_cluster_X.png
    │   ├── subsample_cluster_X.png
    │   ├── Neighborhood_Enrichment_Heatmap_X.png
    │   ├── morphology_image_X.png
    │   └── 2nd_last.png
    ├── Xenium_data_.ipynb
    └── Readme.md
```

---

##  Installation & Setup

```bash
pip install scanpy squidpy anndata matplotlib
```

```bash
jupyter notebook
```

---

##  How to Run

1. Open each module folder
2. Run the notebook step-by-step
3. Refer to module-specific README for explanations

---

##  Future Work

* Integration with machine learning models
* Cell–cell communication analysis
* Multi-omics spatial integration
* 3D spatial transcriptomics

---

## References

* Scanpy Spatial Tutorial
* Squidpy Visium Tutorials
* Xenium Spatial Analysis

---

##  Final Note

This project demonstrates how combining **gene expression + spatial context + imaging** provides deeper biological insight than traditional methods.

---

## Author
Syeda Momina

---
