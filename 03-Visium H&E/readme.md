# MODULE 3: Spatial Transcriptomics Analysis using Squidpy (Visium H&E)

---

## OVERVIEW
This module demonstrates spatial transcriptomics analysis using:
- Squidpy
- Scanpy
- AnnData

Dataset: 10x Genomics Visium (Mouse brain H&E stained section)

We explore:
- Spatial gene expression
- Tissue morphology
- Image-based clustering
- Cell interaction analysis
- Spatial gene variability

---

## LIBRARIES
```python
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import squidpy as sq
```

---

## DATA LOADING

```python
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()
```

### Explanation
- `adata` → gene expression + metadata + spatial coordinates  
- `img` → histology (H&E image)

---

## STEP 1: SPATIAL VISUALIZATION (GENE CLUSTERS)

```python
sq.pl.spatial_scatter(adata, color="cluster")
```

### FIGURE 1: Spatial Cluster Map
<img width="536" height="173" alt="spatial_cluster map H e" src="https://github.com/user-attachments/assets/5cd31dab-1db3-4349-b50b-7c8ec5b48649" />

### Insight
Shows anatomical brain regions based on gene expression (e.g., cortex, hippocampus).

---

## STEP 2: IMAGE FEATURE EXTRACTION

```python
for scale in [1.0, 2.0]:
    sq.im.calculate_image_features(
        adata,
        img.compute(),
        features="summary",
        key_added=f"features_summary_scale{scale}",
        scale=scale,
        n_jobs=4,
    )
```

### Explanation
Extracts morphological features from histology image at multiple scales.

---

## STEP 3: CLUSTERING IMAGE FEATURES

```python
def cluster_features(features):
    adata_tmp = ad.AnnData(features)
    sc.pp.scale(adata_tmp)
    sc.pp.pca(adata_tmp)
    sc.pp.neighbors(adata_tmp)
    sc.tl.leiden(adata_tmp)
    return adata_tmp.obs["leiden"]

adata.obs["features_cluster"] = cluster_features(adata.obsm["features"])
```

---

## STEP 4: COMPARE GENE vs IMAGE CLUSTERS

```python
sq.pl.spatial_scatter(
    adata,
    color=["features_cluster", "cluster"]
)
```

### FIGURE 2: Gene vs Image Clustering
<img width="868" height="244" alt="Gene vs Image Clustering_h e" src="https://github.com/user-attachments/assets/d17b5fae-8f16-4403-b237-cc008a4c58c9" />

### Insight
- Gene-based clusters ≠ always image-based clusters  
- Reveals complementary biological information

---

## STEP 5: SPATIAL CO-OCCURRENCE ANALYSIS

```python
sq.gr.co_occurrence(adata, cluster_key="cluster")

sq.pl.co_occurrence(
    adata,
    cluster_key="cluster",
    clusters="Hippocampus",
)
```

### FIGURE 3: Co-occurrence Heatmap
<img width="458" height="301" alt="Co-occurrence Heatmap_h e" src="https://github.com/user-attachments/assets/5e8aeef9-dd2d-42a9-8730-3f9a8383d6d1" />

### Insight
Shows how frequently clusters appear next to each other spatially.

---

## STEP 6: LIGAND–RECEPTOR INTERACTION ANALYSIS

```python
sq.gr.ligrec(
    adata,
    n_perms=100,
    cluster_key="cluster",
)

sq.pl.ligrec(
    adata,
    cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=["Pyramidal_layer", "Pyramidal_layer_dentate_gyrus"],
)
```

### FIGURE 4: Ligand-Receptor Interactions
*(Place output plot here)*

### Insight
Reveals possible cell communication pathways between brain subregions.

---

## STEP 7: SPATIALLY VARIABLE GENES (MORAN’S I)

```python
genes = adata[:, adata.var.highly_variable].var_names[:1000]

sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
)
```
<img width="425" height="140" alt="moron" src="https://github.com/user-attachments/assets/07765660-2535-4aee-b3e9-81887ec8b46b" />

---

## STEP 8: VISUALIZE SPATIALLY VARIABLE GENES

```python
sq.pl.spatial_scatter(
    adata,
    color=["Olfm1", "Plp1", "Itpka", "cluster"]
)
```

### FIGURE 5: Spatial Gene Expression Maps

<img width="890" height="166" alt="spatial_cluster h e" src="https://github.com/user-attachments/assets/633b5049-fb61-473a-9051-43fd629b4bb7" />


### Insight
Genes show spatially structured expression patterns → indicates non-random biological organization.

---

## FINAL SUMMARY

This module integrates:
- Gene expression (transcriptomics)
- Tissue morphology (H&E imaging)
- Spatial statistics
- Cell-cell communication

### KEY TAKEAWAY
Spatial transcriptomics reveals how gene activity is organized within tissue architecture — something not possible with standard single-cell RNA sequencing.

---

## FIGURE PLACEMENT MAP
- FIGURE 1 → Spatial cluster map  
- FIGURE 2 → Gene vs image clustering  
- FIGURE 3 → Co-occurrence heatmap  
- FIGURE 4 → Ligand-receptor interactions  
- FIGURE 5 → Spatial gene expression maps  

---
