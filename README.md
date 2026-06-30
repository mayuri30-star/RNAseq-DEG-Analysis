# RNAseq-DEG-Analysis
RNA-seq differential gene expression analysis pipeline using DESeq2, edgeR, EnrichR, and DAVID for GEO datasets.

# RNA-seq Differential Gene Expression Analysis

## Overview

This repository performs RNA-seq differential expression analysis using **DESeq2** and **edgeR** on RNA sequencing count data downloaded from the NCBI GEO database.

The workflow identifies differentially expressed genes between two experimental conditions and performs Gene Ontology enrichment analysis using EnrichR and DAVID.

---

## Features

- Build gene count matrix from GEO Excel files
- Differential expression analysis using DESeq2
- Differential expression analysis using edgeR
- PCA visualization
- Volcano plots
- Identification of common significant genes
- GO enrichment using EnrichR
- Functional enrichment using DAVID

---

## Workflow

```
Raw GEO Excel Files
        │
        ▼
Build Count Matrix
        │
        ▼
Sample Metadata
        │
        ▼
DESeq2 Analysis
        │
        ▼
edgeR Analysis
        │
        ▼
Common Differentially Expressed Genes
        │
        ▼
GO Enrichment (EnrichR)
        │
        ▼
DAVID Functional Enrichment
```

---

## Input

The script expects all GEO Excel files to be present inside

```
GSE309456_RAW/
```

Each Excel file should contain:

- Name
- Total exon reads

---

## Software Requirements

R (>=4.3)

Required packages

- ggplot2
- dplyr
- readxl
- DESeq2
- edgeR
- enrichR

Install packages

```r
install.packages(c(
"ggplot2",
"dplyr",
"readxl",
"enrichR"
))

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

BiocManager::install(c(
"DESeq2",
"edgeR"
))
```

---

## Output Files

The analysis generates

### DESeq2

- DESeq2_all_genes.csv
- DESeq2_significant.csv
- DESeq2_UP.csv
- DESeq2_DOWN.csv

### edgeR

- edgeR_all_genes.csv
- edgeR_significant.csv
- edgeR_UP.csv
- edgeR_DOWN.csv

### Figures

- PCA_plot.png
- DESeq2_Volcano.png
- edgeR_Volcano.png

### GO Enrichment

- GO_UP_Top10.csv
- GO_DOWN_Top10.csv

### DAVID

- DAVID_common_UP.txt
- DAVID_common_DOWN.txt

---

## Differential Expression Criteria

Genes were considered significantly differentially expressed using

- Absolute Fold Change > 1.5
- False Discovery Rate (FDR) < 0.1

---

## Analysis Methods

- DESeq2
- edgeR
- EnrichR
- DAVID Functional Annotation
---

## Author

Mayuri Dhakane

M.Sc. Bioinformatics

Savitribai Phule Pune University
