# Single-Cell RNA-seq Analysis of Pancreatic Ductal Adenocarcinoma (PDAC)

## Project Overview

This repository contains an end-to-end single-cell RNA sequencing (scRNA-seq) analysis workflow for **Pancreatic Ductal Adenocarcinoma (PDAC)** using the public **GSE212966** dataset. The project demonstrates a complete analysis pipeline, from raw sequencing data processing to downstream biological interpretation, with the objective of identifying transcriptional programs and genes enriched in malignant epithelial cells while characterizing the tumor microenvironment.

In addition to reproducing a standard single-cell workflow, this project evaluates multiple computational strategies for data integration, cell-type annotation, differential expression analysis, and functional state scoring to improve the robustness of biological interpretation.

## Objectives

* Process raw scRNA-seq data into high-quality gene expression matrices.
* Perform quality control, filtering, and normalization.
* Compare multiple data integration methods to reduce batch effects.
* Identify and annotate major cell populations within the PDAC tumor microenvironment.
* Evaluate and compare automated and marker-based cell annotation approaches.
* Identify genes specifically enriched in malignant epithelial (cancer) cells.
* Perform differential expression analysis at both the single-cell and pseudobulk levels.
* Assess pathway activity using functional state scoring with MSigDB gene sets.

## Dataset

* **Disease:** Pancreatic Ductal Adenocarcinoma (PDAC)
* **Source:** GEO
* **Accession:** GSE212966

## Workflow

1. Raw FASTQ processing using the Galaxy platform
2. Quality control and cell filtering
3. Data normalization and feature selection
4. Batch integration
5. Dimensionality reduction (PCA and UMAP)
6. Cell clustering
7. Cell annotation using:

   * Canonical marker genes
   * SingleR (HPCA)
   * SingleR (BlueprintEncode)
   * scType
   * scCATCH
8. Identification of malignant epithelial cells
9. Cell-level differential expression using Seurat FindMarkers
10. Pseudobulk differential expression analysis
11. Functional state scoring using MSigDB gene sets (including EMT-related pathways)

## Key Results

* Generated integrated single-cell datasets with reduced batch effects.
* Identified distinct immune, stromal, and epithelial cell populations.
* Compared multiple annotation strategies to improve confidence in cell identity assignment.
* Identified genes enriched in malignant epithelial cells.
* Validated differential expression using both cell-level and pseudobulk approaches.
* Characterized pathway activity through single-cell functional state scoring.

## Tools and Packages

* R
* Seurat
* SingleR
* scType
* scCATCH
* Galaxy
* ggplot2
* ComplexHeatmap
* MSigDB gene sets

## Repository Contents

* Data preprocessing scripts
* Quality control workflow
* Integration and clustering analyses
* Cell annotation workflows
* Differential expression analyses
* Functional state scoring
* Figures and visualizations
* Documentation

## Skills Demonstrated

* Single-cell transcriptomics
* Cancer genomics
* Bioinformatics workflow development
* Cell-type annotation
* Differential expression analysis
* Pseudobulk analysis
* Functional pathway scoring
* Data visualization
* Reproducible research
