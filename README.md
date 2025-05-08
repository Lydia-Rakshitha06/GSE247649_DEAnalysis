# GSE247649_DEAnalysis
Differential Gene Expression Analysis using DESeq2 on GSE247649 dataset (Rv vs Ra)

This repository contains R code and output files for a differential gene expression analysis (DEG) of the GSE247649 dataset using the DESeq2 package in R.

## 📌 Project Description

- **Dataset**: GSE247649
- **Comparison**: Rv vs Ra samples
- **Methodology**: DESeq2 pipeline for normalization and DEG analysis, followed by heatmap visualization of significant genes.

## 📂 Files in This Repo

| File | Description |
|------|-------------|
| `GSE247649_RV_RA.R` | Main R script to perform the DEG analysis |
| `GSE247649_Rv.txt` | Raw count data for condition Rv |
| `GSE247649_Ra.txt` | Raw count data for condition Ra |
| `GSE247649_differential_genes.csv` | Output table of significant genes (p < 0.05) |
| `GSE247649_heatmap.png` | Heatmap showing expression patterns of DEGs |
| `README.md` | Project overview and instructions |

## 🛠️ Tools Used

- **R** (v4.0+)
- **DESeq2**
- **ComplexHeatmap**
- **org.Hs.eg.db**
- **ggplot2**, **Cairo**

## 📊 Pipeline Overview

1. Load and merge raw count data
2. Prepare sample metadata
3. Run DESeq2 for normalization and testing
4. Extract significant genes (p < 0.05)
5. Plot heatmap of z-score-normalized DEGs

## 🔧 How to Run

```r
# Ensure all required libraries are installed
source("DESeq2_analysis_GSE247649.R")

