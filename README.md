# Differential Expression Gene (DEG) Analysis - GDS5211

This repository contains an R workflow for analyzing a microarray dataset from NCBI GEO (GDS5211) using the limma pipeline and multiple plotting outputs.

## Overview

The workflow includes:
- Downloading the dataset from NCBI GEO using GEOquery
- Reading and exploring metadata/sample information
- Filtering and normalization of expression data
- Differential expression analysis
- Volcano plot, heatmap, boxplot, PCA, density, and MA plot
- Export of DEG tables and plots

## Project structure

```text
.
├── README.md
├── requirements.R
├── 01_download_and_explore.R
├── 02_preprocessing_normalization.R
├── 03_deg_analysis.R
├── 04_visualization.R
├── 05_full_pipeline.R
├── data/
│   └── .gitkeep
├── results/
│   ├── plots/
│   ├── tables/
│   └── data/
└── .gitignore
```

## Dataset

- Accession: GDS5211
- Source: NCBI GEO
- Platform: microarray

## Running the full pipeline

```r
source("05_full_pipeline.R")
```

## Installing dependencies

```r
source("requirements.R")
```

## Output

The script will generate:
- DEG results table: `results/tables/deg_results_gds5211.csv`
- Volcano plot: `results/plots/volcano_plot.png`
- Heatmap: `results/plots/heatmap_top_genes.png`
- Boxplot: `results/plots/boxplot_top_gene.png`
- PCA plot: `results/plots/pca_plot.png`
- Density plot: `results/plots/density_plot.png`
- MA plot: `results/plots/ma_plot.png`
