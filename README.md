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
Spatial-Transcriptomics-Complete-Pipeline/
│
├── README.md
├── assets/
├── 01_Scanpy_Spatial_Basic/
├── 02_Squidpy_Visium_Fluorescence/
├── 03_Squidpy_Visium_HnE/
├── 04_Squidpy_Xenium/
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
