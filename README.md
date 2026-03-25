# bbb-scRNAseq-visualization
Visualization and summary scripts for single-cell RNA-seq (scRNA-seq) data from mouse brain following BBB-shuttle–mediated siRNA delivery.

This repository reproduces **Figure 3 and associated supplementary figures** from the manuscript.

---
### Overview ###
This repository provides a self-contained R pipeline to:
* Visualize UMAP embeddings by cell type and treatment
* Generate cell-type marker heatmap
* Summarize cell composition across treatment groups
* Plot gene-of-interest expression across cell types (per-cell and pseudobulk)
* Export publication-ready figures and summary tables
---

### Repository Structure ###
```
bbb-scRNAseq-visualization/
├── Figure3_scRNAseq_visualization_pipeline.R
├── Per Cell Data.zip
├── README.md
├── LICENSE
└── .gitignore
```
---

### Data Requirements ###
This repository includes:
* `Per Cell Data.zip` → contains required CSV inputs for plotting
**Important:**
The Seurat object is **not included** due to file size (~1 GB).
---

### Setup Instructions ###
# 1. Download repository
Download this repository from GitHub.
---

# 2. Extract required data
Unzip:
```
Per Cell Data.zip
```
After extraction, your directory should look like:
```
project_dir/
├── Figure3_scRNAseq_visualization_pipeline.R
├── Per Cell Data/
│   ├── DESeq2 comparisions_Hprt - per Cell Type.csv
│   ├── Pseudobulked expression data_Hprt - per Cell Type.csv
│   └── Pseudobulk Normalized Counts_All Genes - per Cell Type.csv
```
---
# 3. Provide Seurat object
Place the processed Seurat object in your project directory:
```
SCTint_srat_25-1b23_QC_filtered.rds
```
This file will be available via:
* GEO (accession to be provided upon publication)
---

# 4. Set working directory
In the script:
```r
project_dir <- "path/to/your/project_folder"
```
---
# 5. Run the pipeline
---

### Outputs ###

All outputs are saved to:
```
outputs/
├── figures/
│   ├── UMAP/
│   ├── Markers/
│   ├── Barplots/
│   └── Pseudobulk/
├── tables/
│   ├── Summary/
│   └── Pseudobulk/
```

Generated figures include:

* UMAP (cell type + treatment)
* Marker heatmap
* Cell composition summaries
* Gene expression barplots (free + locked scaling)
* Pseudobulk gene panels
---

### Optional Customization ###
Change gene of interest
```r
target_gene <- "Tfrc"
```
Replace with any gene present in the dataset.

---

 
## R packages required:

* Seurat
* dplyr
* ggplot2
* tidyr
* stringr
* ggpubr
* purrr
* readr
* scales
* grid
* tibble
* Matrix


### Data Availability ###
* Raw and processed scRNA-seq data: **GEO (accession pending)**
* Large objects (Seurat): available via GEO or upon request
---

### License ###
MIT License
---

Prepared by:
**Lynette Bustos**

---
###  Notes ###
* This repository contains **visualization-focused code only**
* Upstream processing (QC, clustering, annotation) was performed separately
* Scripts are designed for **reproducibility of manuscript figures**
---
