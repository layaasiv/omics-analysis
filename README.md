# Single-Cell & Spatial Omics Analysis Workflows

End-to-end analysis workflows for single-cell RNA-seq, single-cell ATAC-seq, and spatial transcriptomics data, implemented in R.

## Tools & Technologies
- **scRNA-seq:** Seurat
- **scATAC-seq:** Signac
- **Spatial transcriptomics:** Seurat (Visium)
- **Language:** R
- **Data sources:** GEO, 10x Genomics

## Overview
This repository contains independent analysis workflows for three single-cell and spatial omics modalities. Each workflow follows best practices for its respective data type, from raw data ingestion through dimensionality reduction, clustering, and visualization.

## Workflows

### scRNA-seq | Seurat
Analysis of scRNA-seq data from hepatoblastoma samples (3 individuals).
- QC filtering (nFeature, nCount, mitochondrial percentage thresholds)
- Normalization and highly variable gene selection
- Dimensionality reduction (PCA, UMAP)
- Graph-based clustering
- Batch correction and integration across 3 donors
- Marker gene identification and cluster annotation
- Visualization: UMAP embeddings, feature plots, dot plots, violin plots

### scATAC-seq | Signac
Analysis of single-cell chromatin accessibility data.
- Fragment file processing and peak quantification
- QC filtering using TSS enrichment score and nucleosome signal
- Dimensionality reduction (LSI, UMAP)
- Graph-based clustering
- Integration with matched scRNA-seq data
- Chromatin accessibility visualization: coverage plots, genomic region tracks

### Spatial Transcriptomics | Seurat (Visium)
Analysis of 10x Visium spatial transcriptomics data.
- Spatial data loading and preprocessing
- Normalization and dimensionality reduction
- Spatially variable gene identification
- Visualization of gene expression in tissue context

## Setup
```r
install.packages("Seurat")
# For Signac:
install.packages("Signac")
# Required Bioconductor packages:
BiocManager::install(c("GenomicRanges", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38"))
```
