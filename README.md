# MiniProject_Microarray_RMA_Normalization_GSE148537
R based microarray pipeline on GSE148537 breast cancer data using GEOquery and affy. Performed RMA normalization with background correction, quantile normalization, and log2 summarization, generated expression matrix, mapped probe IDs to gene symbols, and enabled downstream differential expression analysis

---

# Breast Cancer Microarray Preprocessing and Annotation Pipeline (GSE148537)

## Overview

This mini project implements a reproducible R based microarray preprocessing workflow on the human breast cancer dataset **GSE148537** obtained from the **Gene Expression Omnibus**. The objective was to process raw CEL files, perform robust normalization, and generate a biologically annotated gene expression matrix suitable for downstream transcriptomic analysis.

## Technical Workflow

### 1. Data Acquisition

Supplementary RAW CEL files were downloaded using `getGEOSuppFiles()` from GEOquery and extracted from the compressed TAR archive.

### 2. RMA Normalization

Raw microarray intensity data were imported using `ReadAffy()` from the affy package. Robust Multi array Average (RMA) normalization was applied, which includes:

* Background correction
* Quantile normalization
* Log2 transformation and probe level summarization

This generated a normalized expression matrix with log2 scaled values typically ranging from 2 to 14.

### 3. Probe Annotation

Processed expression data were converted into a structured dataframe. Probe IDs were mapped to gene symbols by extracting feature metadata from the GEO ExpressionSet object using `fData()`. Annotation merging was performed using `tibble::rownames_to_column()` and `dplyr::inner_join()` to generate a biologically interpretable gene level matrix.

## Output

The final output is an annotated gene expression dataframe linking probe IDs to official gene symbols, enabling downstream analyses such as differential expression testing, log2 fold change estimation, clustering, pathway enrichment, and candidate biomarker discovery.

## Tools and Packages

* affy
* GEOquery
* tidyverse
* Bioconductor

This workflow ensures reproducible preprocessing of microarray data and provides a clean foundation for advanced breast cancer transcriptomic analysis.

