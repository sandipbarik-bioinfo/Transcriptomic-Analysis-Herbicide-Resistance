# Transcriptomic Analysis of Pinoxaden Resistance in *Apera spica-venti*

RNA-seq transcriptomic analysis investigating gene-expression changes associated with the development of pinoxaden resistance in *Apera spica-venti*.

## Overview

This repository contains the R-based transcriptomic analysis workflow used to investigate differential gene expression associated with herbicide resistance in *Apera spica-venti*.

The analysis compares transcriptomic profiles between herbicide-resistant and susceptible plants under control and herbicide-treated conditions. The workflow includes differential expression analysis, identification of significantly differentially expressed genes (DEGs), comparison of expression patterns between experimental conditions, gene-level analysis, and visualization of transcriptomic patterns.

The project was developed to investigate molecular responses associated with herbicide resistance and to identify genes and expression patterns that may contribute to the resistant phenotype.

---

## Project Objectives

The main objectives of this project were to:

- Analyze transcriptomic data from resistant and susceptible plant populations.
- Compare gene expression between resistant and susceptible conditions.
- Identify significantly differentially expressed genes (DEGs).
- Apply statistical and fold-change based filtering to prioritize biologically relevant genes.
- Compare DEG sets between experimental conditions.
- Identify common and condition-specific genes.
- Characterize genes using gene identifiers and functional annotation resources.
- Visualize global transcriptomic patterns using PCA, Venn diagrams, and other graphical analyses.
- Generate gene-level datasets for downstream biological interpretation.

---
## Dataset

The primary expression dataset used in this analysis was obtained from the NCBI Gene Expression Omnibus (GEO), accession **GSE204857**.

The dataset contains transcriptomic data from susceptible and resistant *Apera spica-venti* populations under control and pinoxaden-treated conditions.

---

## Differential Expression Analysis

Differential expression analysis was performed using an R-based workflow.

Multiple experimental comparisons were analyzed to investigate transcriptional differences between resistant and susceptible plants under control and herbicide-treated conditions.

Genes were considered significantly differentially expressed based on:

- Adjusted *p*-value < 0.05
- Absolute log2 fold change ≥ 1

These criteria were used to focus the analysis on genes showing statistically significant and biologically relevant changes in expression.

---

## Analysis Workflow

The computational workflow consisted of the following major stages:

1. Transcriptomic data preparation and preprocessing
2. Quality assessment and exploratory analysis
3. Differential expression analysis
4. Statistical and fold-change filtering of DEGs
5. Comparison of DEG sets between experimental conditions
6. Identification of common and condition-specific genes
7. Gene identifier and gene-symbol analysis
8. Functional/gene annotation
9. Visualization of transcriptomic patterns
10. Generation of final gene-level datasets for biological interpretation
