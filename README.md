# MEANDER
### Revealing the mechanisms of aneuploidy-driven cell fate decisions in early human development

One of the MEANDER work packages develops a framework for **reference-independent copy number variation (CNV) inference from single-cell RNA-sequencing data**, being part of a broader research initiative investigating **aneuploidy-driven cellular dynamics during early human development**. 

---

# Overview

Chromosomal mosaicism is common in developmental systems. However, most existing **single-cell RNA-seq CNV inference tools rely on predefined diploid reference populations** to detect copy number alterations from gene expression patterns.

This requirement presents challenges in systems where:

- a **true diploid reference population may not exist**
- the **reference population may itself be mosaic**
- experimental designs **lack clear control populations**

Several reference-free approaches have been proposed, but these methods often rely on **stringent thresholds** and typically detect **whole-chromosome events**, limiting their ability to identify **segmental CNVs** and reducing resolution.

**MEANDER** aims to address these limitations by introducing a **cross-fold reference strategy** for CNV inference.

# Input Data

The current preliminary MEANDER pipeline on scRNA-seq datasets containing:
- gene expression matrix
- gene genomic coordinates
- cell metadata

---


# Implementation

The underdevelopment framework is implemented in **R** and designed as a modular pipeline for CNV inference from single-cell RNA sequencing datasets.

## Required R Packages

### CRAN packages

```r
cran_packages <- c(
  "dplyr",
  "tidyr",
  "data.table",
  "ggplot2",
  "patchwork",
  "cowplot",
  "purrr",
  "igraph",
  "Seurat",
  "msigdbr"
)

bioc_packages <- c(
  "edgeR",
  "GSVA",
  "GenomicRanges",
  "IRanges"
)

install.packages(cran_packages)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(bioc_packages)
```

---

# Repository Status

⚠️ **Active Development**

This repository is maintained as part of an ongoing research project and is currently under active development.

Some components may:

- be under refinement
- contain experimental code
- be reorganized as the framework matures

The repository is shared to ensure **transparency**


# Authors

Pedro Granjo and  Claudia Spits from the VUB - Genetics Reproduction & Development Research Group

---

# Contact

For questions regarding the project or collaboration opportunities:
Claudia Spits  
Email: Claudia.Spits@vub.be