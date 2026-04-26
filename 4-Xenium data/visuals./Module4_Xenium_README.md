# 🧬 Module 4: Spatial Transcriptomics with Xenium In Situ Data (Squidpy + SpatialData)

---

## 🎯 Objective

This module demonstrates the analysis of **10x Genomics Xenium in situ transcriptomics data** using Squidpy integrated with the SpatialData framework.

Unlike Visium (spot-based, ~50µm resolution), Xenium provides **true single-cell resolution** gene expression measured directly within intact tissue — no deconvolution, no ambiguity. This dataset integrates:

* Single-cell gene expression profiles (480 genes across 161,000 cells)
* High-resolution morphology focus image of lung tissue
* Cell and nucleus boundary segmentations
* Individual transcript coordinates (~40 million points in 3D)

This allows us to study **how individual cells are spatially organized within human lung cancer tissue** and uncover cell-cell interaction patterns invisible to bulk or spot-based methods.

---

## 📦 Dataset Overview

We use the **10x Genomics Xenium Human Lung Cancer dataset** (`Xenium_V1_human_Lung_2fov`), downloaded directly from 10x Genomics.

### SpatialData Object Structure:

| Layer | Contents |
|---|---|
| **Images** | `morphology_focus` — multi-scale grayscale tissue image |
| **Labels** | `cell_labels`, `nucleus_labels` — pixel-level segmentation masks |
| **Points** | `transcripts` — 40M+ individual transcript detections (x, y, z) |
| **Shapes** | `cell_boundaries`, `nucleus_boundaries`, `cell_circles` |
| **Tables** | `table` — AnnData (161,000 cells × 480 genes) |

📌 Every element shares a common coordinate system (`global`), making it possible to overlay expression, segmentation, and morphology images precisely.

---

## 🔄 Workflow Overview

```
Download Dataset (wget → unzip)
        ↓
Load with spatialdata-io (xenium reader)
        ↓
Quality Control (control probe % + 4-panel histograms)
        ↓
Filter Cells & Genes
        ↓
Normalize → Log → PCA → Neighbors → UMAP → Leiden Clustering
        ↓
Visualize: UMAP + Spatial Scatter
        ↓
Build Spatial Neighborhood Graph (Delaunay)
        ↓
Centrality Scores
        ↓
Co-occurrence Probability
        ↓
Neighborhood Enrichment
        ↓
Moran's I Spatial Autocorrelation
        ↓
Gene Expression on Morphology Image (spatialdata-plot)
        ↓
Interactive Visualization (napari-spatialdata)
```

---

# ⚙️ STEP 0 — Setup & Installation

```python
!pip install spatialdata spatialdata-io
!pip install squidpy scanpy seaborn matplotlib
!pip install leidenalg igraph
```

```python
import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
```

> 📌 `leidenalg` and `igraph` are required for Leiden clustering. `spatialdata-io` provides the Xenium reader. Install all packages before running the notebook.

---

# 📥 STEP 1 — Download & Load Data

The Xenium dataset is downloaded directly from 10x Genomics and extracted into the working directory:

```python
!wget https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Lung_2fov/Xenium_V1_human_Lung_2fov_outs.zip

import zipfile
with zipfile.ZipFile("Xenium_V1_human_Lung_2fov_outs.zip", "r") as z:
    z.extractall("/content")
```

We then load all data layers using the `spatialdata-io` Xenium reader:

```python
xenium_path = "/content"
sdata = xenium(xenium_path)
```

The AnnData expression table is extracted for all downstream analysis:

```python
adata = sdata.tables["table"]
# AnnData: 161,000 cells × 480 genes
```

Each row in `adata.obs` is one segmented cell. Key columns:

| Column | Meaning |
|---|---|
| `transcript_counts` | Raw Xenium transcript detections per cell |
| `total_counts` | Total counts including controls |
| `cell_area` | Physical area of the segmented cell (µm²) |
| `nucleus_area` | Physical area of the nucleus (µm²) |
| `z_level` | Z-plane of the cell |
| `nucleus_count` | Number of nuclei detected per cell |

Cell spatial coordinates are stored in `adata.obsm["spatial"]` as (x, y) pairs — automatically set by the Xenium reader.

> 📌 **Why SpatialData?** The SpatialData framework stores all modalities (images, segmentation, transcripts, expression table) in a single object with a shared coordinate system. This makes spatial overlays exact and reproducible across all analysis steps.

---

# 📊 STEP 2 — Quality Control

We compute the fraction of control probe and control codeword detections — these are spike-in negative controls that should be near zero in a high-quality Xenium run:

```python
cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")
```

Expected output:
```
Negative DNA probe count % : ~0.005%
Negative decoding count % : ~0.003%
```

Near-zero values confirm **very low background noise** — the Xenium assay is highly specific with minimal false detections.

---

## 📸 VISUAL 1 — QC Distributions (4-Panel Histogram)

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run the block below, save the matplotlib figure, and upload it.
> **Suggested filename:** `qc_distributions_xenium.png`

```python
fig, axs = plt.subplots(1, 4, figsize=(15, 4))

axs[0].set_title("Total counts per cell")
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])

axs[1].set_title("Transcript counts per cell")
sns.histplot(adata.obs["transcript_counts"], kde=False, ax=axs[1])

axs[2].set_title("Cell area")
sns.histplot(adata.obs["cell_area"], kde=False, ax=axs[2])

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)
```

| Panel | Metric | What it tells you |
|---|---|---|
| 1 | **Total counts per cell** | Overall RNA capture — very low values = empty/dead cells |
| 2 | **Transcript counts per cell** | Directly measured Xenium transcripts per cell |
| 3 | **Cell area** | Segmented cell size — outliers may be debris or doublets |
| 4 | **Nucleus ratio** (nucleus_area / cell_area) | Nuclear-to-cytoplasm ratio — varies systematically by cell type |

📌 **Insight:** These four distributions directly inform the filtering thresholds applied in Step 3. Cells with very few counts (far left tail in panels 1–2) are low-quality and should be removed. Cells with unusually large area may be segmentation artifacts.

---

# 🧹 STEP 3 — Filtering & Normalization

Based on the QC plots, we remove low-quality cells and rarely-detected genes:

```python
sc.pp.filter_cells(adata, min_counts=10)   # Remove cells with fewer than 10 total counts
sc.pp.filter_genes(adata, min_cells=5)     # Remove genes detected in fewer than 5 cells
```

We then save raw counts and normalize for downstream analysis:

```python
adata.layers["counts"] = adata.X.copy()   # Preserve raw integers for future tools

sc.pp.normalize_total(adata, inplace=True) # Scale each cell to equal total count sum
sc.pp.log1p(adata)                         # Log-transform to stabilize variance
```

> 📌 **Why save `adata.layers["counts"]`?** The raw count layer is required by differential expression tools, trajectory analysis, and spatial deconvolution methods that expect unnormalized integer counts as input.

---

# 🔬 STEP 4 — Dimensionality Reduction & Leiden Clustering

```python
sc.pp.pca(adata)         # Reduce 480 genes → principal components
sc.pp.neighbors(adata)   # Build k-nearest-neighbor graph in PCA space
sc.tl.umap(adata)        # Non-linear 2D embedding for visualization
sc.tl.leiden(adata)      # Graph-based community detection → cluster labels
```

This standard Scanpy pipeline compresses the 480-dimensional gene space into a 2D UMAP for visualization and identifies transcriptomically distinct cell populations (Leiden clusters) using graph community detection.

---

## 📸 VISUAL 2 — UMAP: Total Counts & Leiden Clusters

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run `sc.pl.umap(...)` and save the output.
> **Suggested filename:** `umap_leiden_xenium.png`

```python
sc.pl.umap(
    adata,
    color=["total_counts", "leiden"],
    wspace=0.4,
)
```

📌 **Insight:**
- The **total counts** panel checks for systematic quality differences between clusters — a cluster composed entirely of low-count cells is likely a technical artifact.
- The **Leiden** panel reveals distinct transcriptomic populations. In a lung cancer dataset, these typically include tumor epithelial cells, T cells, macrophages, fibroblasts, endothelial cells, and B cells — each with a characteristic gene signature.

---

## 📸 VISUAL 3 — Spatial Scatter: Leiden Clusters in Tissue

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run `sq.pl.spatial_scatter(adata, ...)` and save the output.
> **Suggested filename:** `spatial_leiden_xenium.png`

```python
sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=["leiden"],
    wspace=0.4,
)
```

📌 **Insight:** This maps Leiden cluster identity back onto actual tissue X-Y coordinates. Unlike the abstract UMAP, this reveals **where each cell type physically lives** in the tumor section — for example, whether immune cells infiltrate the tumor core or concentrate at the tumor-stroma boundary.

---

# 🕸️ STEP 5 — Build Spatial Neighborhood Graph

All spatial statistics require a connectivity graph that defines which cells are neighbors. We build this from spatial coordinates using Delaunay triangulation:

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```

> 📌 **Why Delaunay?** For single-cell data where cells are irregular in position and density (unlike the regular Visium grid), Delaunay triangulation builds a geometrically natural graph — each cell connects to its true nearest neighbors without needing a fixed search radius. The graph adapts automatically to local tissue density.

The graph is stored in `adata.obsp["spatial_connectivities"]` and used by Steps 6–9.

---

# 📐 STEP 6 — Centrality Scores

Centrality scores describe the **topological role of each cluster** in the spatial neighborhood graph — are cells of this type central hubs of the tissue, or spatially isolated?

```python
sq.gr.centrality_scores(adata, cluster_key="leiden")
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))
```

| Score | Definition | High value means... |
|---|---|---|
| **Closeness centrality** | Average inverse distance to all other nodes | Cluster is at the spatial center of the tissue |
| **Degree centrality** | Fraction of cross-cluster direct connections | Cluster interfaces with many other cell types |
| **Clustering coefficient** | How tightly the cluster's own cells interconnect | Cluster forms a dense, cohesive spatial patch |

---

## 📸 VISUAL 4 — Centrality Scores per Leiden Cluster

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run `sq.pl.centrality_scores(...)` and save the output.
> **Suggested filename:** `centrality_scores_xenium.png`

📌 **Insight:** Clusters with high **closeness centrality** are likely stromal or immune cell types that act as **spatial hubs** — bridging different tissue compartments. Clusters with high **clustering coefficient** form tight, cohesive patches — characteristic of tumor nests, lymphoid follicles, or glandular structures.

---

# 🔁 STEP 7 — Co-occurrence Probability

The co-occurrence score answers: **"Given that cluster X is present at a location, how much more (or less) likely is cluster Y to appear within radius r?"**

$$\text{score}(r) = \frac{p(\text{cluster Y} \mid \text{cluster X at radius } r)}{p(\text{cluster Y})}$$

Score > 1 → spatial attraction at that radius. Score < 1 → spatial exclusion. Computed across a range of increasing radii.

Because co-occurrence is computationally intensive, we first create a 50% subsample:

```python
sdata.tables["subsample"] = sc.pp.subsample(adata, fraction=0.5, copy=True)
adata_subsample = sdata.tables["subsample"]

sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)

sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",       # Focus: co-occurrence around cluster 12
    figsize=(10, 10),
)

sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
)
```

---

## 📸 VISUAL 5 — Co-occurrence Probability: Cluster 12 vs All Clusters

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run `sq.pl.co_occurrence(..., clusters="12")` and save the output.
> **Suggested filename:** `co_occurrence_cluster12_xenium.png`

📌 **Insight:** Each line in the plot is a different Leiden cluster. Lines rising above 1.0 at short radii indicate clusters **preferentially close to cluster 12** — these are its direct spatial interaction partners in the tumor microenvironment.

---

## 📸 VISUAL 6 — Spatial Scatter: Subsampled Leiden Clusters

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run `sq.pl.spatial_scatter(adata_subsample, color="leiden", ...)` and save the output.
> **Suggested filename:** `spatial_scatter_subsample_xenium.png`

📌 **Insight:** This is the spatial reference map for interpreting the co-occurrence plot. By comparing which clusters appear adjacent to cluster 12 on the tissue, you can visually validate the statistical co-occurrence signal — building biological intuition about shared microenvironmental niches.

---

# 🏘️ STEP 8 — Neighborhood Enrichment Analysis

Neighborhood enrichment asks: **"Are two clusters significantly more (or less) likely to be direct spatial neighbors than random chance predicts?"**

It computes observed co-localization counts and compares them against 1,000 permutations of shuffled cluster labels, outputting a z-score for each cluster pair:

```python
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
# 1,000 permutations run to establish the null distribution
```

Results are visualized as a symmetric heatmap alongside the spatial map:

```python
fig, ax = plt.subplots(1, 2, figsize=(13, 7))

sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])
```

---

## 📸 VISUAL 7 — Neighborhood Enrichment Heatmap + Spatial Reference

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run the `fig, ax = plt.subplots(1, 2, ...)` block and save the combined figure.
> **Suggested filename:** `nhood_enrichment_xenium.png`

📌 **Insight:**
- **Warm colors (high positive z-scores)** → clusters that are significantly co-localized — likely **cell-cell interaction partners** (e.g., tumor cells and their associated macrophages or T cells).
- **Cool colors (negative z-scores)** → clusters that spatially exclude each other — they occupy distinct tissue compartments and do not directly interact.
- Side-by-side with the spatial map, you can visually confirm which statistically enriched pairs are anatomically adjacent in the tissue.

---

# 📡 STEP 9 — Moran's I Spatial Autocorrelation

Moran's I measures whether a **gene's expression is spatially structured** or randomly scattered across cells:

| Moran's I | Interpretation |
|---|---|
| Close to **+1** | Expression is spatially clustered — similar cells are near each other |
| Close to **0** | Expression is spatially random |
| **Negative** | Expression is overdispersed — similar cells avoid each other |

```python
sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)

sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)

adata_subsample.uns["moranI"].head(10)
```

### Top Spatially Autocorrelated Genes:

| Gene | Moran's I | Biological Role |
|---|---|---|
| **AREG** | ~0.70 | Amphiregulin — EGFR ligand, drives epithelial tumor growth |
| **MET** | ~0.68 | MET receptor tyrosine kinase — invasion and metastasis oncogene |
| **ANXA1** | ~0.67 | Annexin A1 — anti-inflammatory signal, macrophage marker |
| **EPCAM** | ~0.63 | Epithelial Cell Adhesion Molecule — epithelial/tumor cell marker |
| **DMBT1** | ~0.59 | Deleted in Malignant Brain Tumors 1 — tumor suppressor |

📌 **Genes with high Moran's I are the most spatially informative** — their expression reflects tissue architecture and are strong candidates for spatial biomarkers of the tumor microenvironment.

---

## 📸 VISUAL 8 — Spatial Expression of Top Moran's I Genes (AREG & MET)

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run `sq.pl.spatial_scatter(color=["AREG","MET"])` and save the output.
> **Suggested filename:** `moran_genes_spatial_xenium.png`

```python
sq.pl.spatial_scatter(
    adata_subsample,
    library_id="spatial",
    color=["AREG", "MET"],
    shape=None,
    size=2,
    img=False,
)
```

📌 **Insight:** AREG and MET expression appears in contiguous spatial patches — consistent with their high Moran's I. This spatial clustering reflects **clonal tumor expansion**: groups of genetically similar tumor cells proliferating together maintain transcriptional coherence across tissue space.

---

# 🖼️ STEP 10 — Gene Expression Overlaid on Morphology Image

Using `spatialdata-plot`, we render gene expression directly on top of the raw morphology focus image — anchoring molecular findings in tissue histology:

```python
import spatialdata_plot

gene_name = ["AREG", "MET"]
for name in gene_name:
    sdata.pl.render_images("morphology_focus").pl.render_shapes(
        "cell_circles",
        color=name,
        table_name="table",
        use_raw=False,
    ).pl.show(
        title=f"{name} expression over Morphology image",
        coordinate_systems="global",
        figsize=(10, 5),
    )
```

---

## 📸 VISUAL 9 — AREG Expression on Morphology Image

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run the `spatialdata_plot` loop, save the AREG figure.
> **Suggested filename:** `areg_morphology_xenium.png`

📌 **Insight:** Overlaying AREG on the morphology image reveals **which anatomical structures express this tumor growth factor** — for example, glandular epithelial regions versus stromal areas. This is the defining advantage of Xenium: every expression value is physically anchored to a cell whose morphological context is directly visible.

---

## 📸 VISUAL 10 — MET Expression on Morphology Image

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Run the `spatialdata_plot` loop, save the MET figure.
> **Suggested filename:** `met_morphology_xenium.png`

📌 **Insight:** MET focal high-expression regions on the morphology image can reveal areas of oncogenic amplification within the tumor. MET overexpression is a known driver of invasion and metastasis in lung cancer — its spatial pattern identifies which specific cells and tissue regions are most affected.

---

# 🎛️ STEP 11 — Interactive Visualization (napari-spatialdata)

For exploratory analysis, all SpatialData layers are browsable interactively in napari:

```python
from napari_spatialdata import Interactive
Interactive(sdata)
```

napari-spatialdata lets you:
- Pan and zoom across the full tissue section (50,000+ pixels wide)
- Toggle between gene expression, cluster labels, segmentation masks, and raw morphology
- Overlay individual transcript points, color-coded by gene identity
- Interactively select and inspect individual segmented cells

---

## 📸 VISUAL 11 — napari Interactive Viewer (AREG Expression)

> 🖼️ **UPLOAD YOUR IMAGE HERE**
> Take a screenshot of the napari viewer with AREG expression displayed across the tissue.
> **Suggested filename:** `napari_areg_xenium.png`

📌 **Insight:** napari provides a "Google Maps"-style exploration experience for spatial omics data. The ability to zoom from a whole-tissue overview down to individual cell segmentation boundaries — with expression and transcript points overlaid — is essential for hypothesis generation and spatial validation at scale.

---

# 🗂️ Complete Image Upload Checklist

| # | Suggested Filename | Code Block to Run | Step |
|---|---|---|---|
| 1 | `qc_distributions_xenium.png` | `fig, axs = plt.subplots(1, 4, ...)` QC histogram block | Step 2 |
| 2 | `umap_leiden_xenium.png` | `sc.pl.umap(color=["total_counts","leiden"])` | Step 4 |
| 3 | `spatial_leiden_xenium.png` | `sq.pl.spatial_scatter(adata, color=["leiden"])` | Step 4 |
| 4 | `centrality_scores_xenium.png` | `sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16,5))` | Step 6 |
| 5 | `co_occurrence_cluster12_xenium.png` | `sq.pl.co_occurrence(..., clusters="12", figsize=(10,10))` | Step 7 |
| 6 | `spatial_scatter_subsample_xenium.png` | `sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2)` | Step 7 |
| 7 | `nhood_enrichment_xenium.png` | `fig, ax = plt.subplots(1, 2, figsize=(13, 7))` combined block | Step 8 |
| 8 | `moran_genes_spatial_xenium.png` | `sq.pl.spatial_scatter(color=["AREG","MET"], img=False)` | Step 9 |
| 9 | `areg_morphology_xenium.png` | `spatialdata_plot` render loop — AREG output | Step 10 |
| 10 | `met_morphology_xenium.png` | `spatialdata_plot` render loop — MET output | Step 10 |
| 11 | `napari_areg_xenium.png` | Screenshot of napari viewer | Step 11 |

---

# 💡 KEY INSIGHTS

### 🔹 1. Single-cell resolution changes everything
Xenium assigns transcripts to individual cells — no spot mixing, no deconvolution uncertainty. Every data point is a real, segmented, spatially located cell with a known position in tissue.

### 🔹 2. SpatialData scales to real tissue complexity
With 161,000 cells, 40M+ transcripts, and images up to 51,000 pixels wide, the Zarr-backed SpatialData framework is essential — traditional in-memory approaches would fail at this scale.

### 🔹 3. Spatial statistics reveal tumor microenvironment organization
Centrality scores, co-occurrence probabilities, and neighborhood enrichment together reveal **which cell types interact** and **where those interactions occur** — the functional geography of the tumor.

### 🔹 4. Moran's I identifies spatially structured genes
AREG and MET are not just highly expressed — they are **spatially clustered**. Their expression patterns are shaped by tissue architecture, not random cell-to-cell noise. These are the most biologically meaningful spatial biomarkers.

### 🔹 5. Morphology + expression = full biological picture
Overlaying gene expression on the morphology focus image bridges transcriptomics with classical histopathology — connecting molecular discoveries to the tissue structures a pathologist would recognize.

---

# 📁 Final Output Summary

| Object | Contents |
|---|---|
| `sdata` | Full SpatialData object — images, labels, shapes, points, tables |
| `adata` | Filtered, normalized, clustered AnnData for all cells |
| `adata_subsample` | 50% subsample used for co-occurrence and Moran's I |
| `adata.uns["moranI"]` | Moran's I scores and p-values for all 480 genes |
| `adata.obs["leiden"]` | Cluster label for each cell |
| `adata.obsm["X_umap"]` | 2D UMAP embedding coordinates |
| `adata.obsp["spatial_connectivities"]` | Delaunay spatial neighborhood graph |
