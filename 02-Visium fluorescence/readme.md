# 🧬 Module 2: Spatial Transcriptomics with Fluorescence Imaging (Squidpy)

---

##  Objective

This module demonstrates the analysis of **Visium spatial transcriptomics data combined with fluorescence imaging** using Squidpy.

Unlike basic Scanpy workflows, this dataset integrates:

* Gene expression profiles
* High-resolution fluorescence microscopy images
* Image-derived morphological features

This allows us to study how **cellular morphology and gene expression jointly define tissue structure**.

---

## Dataset Overview

We use a **10x Genomics Visium mouse brain dataset**, which includes:

* Spatial gene expression matrix (AnnData)
* Fluorescence microscopy image (ImageContainer)
* Precomputed cluster annotations

###  Fluorescence Channels:

* DAPI → nuclei (DNA)
* NEUN → neurons
* GFAP → glial cells

📌 This multi-channel structure enables **cell-type–specific spatial interpretation**.

---

##  Workflow Overview

```id="flow2"
Load Data → Visualize Tissue → Image Processing → Segmentation → Feature Extraction → Clustering → Spatial Interpretation
```

---

#  STEP 1 — Load Data & Initial Visualization

We begin by loading:

* AnnData (gene expression)
* ImageContainer (fluorescence image)

---

##  VISUAL 1 — Spatial Cluster Overview

📌 Function:

```python id="sq1"
sq.pl.spatial_scatter(adata, color="cluster")
```
<img width="354" height="260" alt="spatial_cluster_flu" src="https://github.com/user-attachments/assets/ed6f5de8-618c-462c-b930-63943da52375" />


📌 Insight:
Shows initial anatomical structure of brain regions in spatial coordinates.

---

##  VISUAL 2 — Fluorescence Channels
📌 Function:

```python id="sq2"
img.show(channelwise=True)
```

<img width="532" height="188" alt="Fluorescence Channels " src="https://github.com/user-attachments/assets/ff7d66f3-4b3c-4294-9d49-43bb33d3320d" />

📌 Insight:
Each channel highlights a different biological structure:

* DAPI → nuclei distribution
* NEUN → neuronal density
* GFAP → glial localization

👉 This connects **molecular biology + imaging**

---

# Image Preprocessing

We apply:

* smoothing (noise reduction)
* segmentation (cell detection)

```python id="sq3"
sq.im.process(img, layer="image", method="smooth")
sq.im.segment(img, method="watershed", channel=0)
```

---

## Cell Segmentation Result

📌 Function:

```python id="sq4"
img.show(layer="segmented_watershed")
```

<img width="491" height="248" alt="img_seg_flu" src="https://github.com/user-attachments/assets/b5769b16-e467-4200-9488-5844ea776555" />


📌 Insight:
Each segmented object corresponds to a **cell nucleus or region**

---

#  Image Feature Extraction 
We then extract quantitative features from tissue images:

* summary (intensity stats)
* histogram (pixel distribution)
* texture (spatial variation)
* segmentation (cell structure)

---

##  Feature Types Comparison

📌 Function:

```python id="sq5"
sq.pl.spatial_scatter(... features_summary_cluster ...)
sq.pl.spatial_scatter(... features_histogram_cluster ...)
sq.pl.spatial_scatter(... features_texture_cluster ...)
```


📌 Insight:
Different image features capture different biological signals:

* intensity → expression richness
* texture → tissue structure
* segmentation → cell density

---

#  Cell-Level Biological Interpretation

We now link:

* image features
* gene expression
* cluster labels

  
<img width="419" height="338" alt="Feature Types Comparison_flu" src="https://github.com/user-attachments/assets/bf6a0b90-9cac-4a01-94cc-91f7a5ad3c44" />

---



# 📊 Multi-scale Feature Clustering 

We combine multiple feature scales:

* original resolution
* contextual crops
* low-resolution views

Then cluster them using Leiden algorithm.

---

##  Feature-Based Clusters vs Gene Clusters

📌 Function:

```python id="sq7"
sq.pl.spatial_scatter(
    adata,
    color=[
        "features_summary_cluster",
        "features_histogram_cluster",
        "features_texture_cluster",
        "cluster"
    ]
)
```

<img width="661" height="316" alt="features vs clusters_flu" src="https://github.com/user-attachments/assets/275c2939-86c2-496f-9cdc-62582875cd79" />


📌 Insight:

* Image-based clusters reveal **substructure not visible in gene-only clustering**
* Hippocampus splits into finer regions
* Cortex shows layered organization

---

#  KEY INSIGHTS

This module demonstrates that:

### 🔹 1. Gene expression alone is not enough

Spatial + image features add missing biological context

### 🔹 2. Fluorescence imaging improves interpretation

Cell types become visually and spatially identifiable

### 🔹 3. Image features reveal hidden structure

Even within the same gene cluster

---

# 📁 FINAL OUTPUT FILES


---

#  FINAL TAKEAWAY

This workflow shows a key concept in modern spatial omics:

>  “Gene expression + imaging + spatial context = true tissue biology”

---
