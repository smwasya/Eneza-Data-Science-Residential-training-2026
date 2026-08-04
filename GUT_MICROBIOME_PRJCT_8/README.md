# Gut Microbiome Project (Project 8)

This repository contains all work related to **Project 8** of the **ENEZA Data Science Initiative (DSI) Residential Training**. 

The project focuses on reproducible 16S rRNA gut microbiome analysis and the development of machine learning models for disease classification using microbial community profiles.

## Repository Structure

```text
GUT_MICROBIOME_PRJCT_8/
├── Data/
├── ENEZAflow/
├── Machine-Learning/
├── Results/
└── Scripts/
```

## Directory Structure

### `Data/`

Contains processed datasets used throughout the project, including:

* Diversity analysis tables
* ASVs tables for downstream analyses
* Input datasets for machine learning models

### `ENEZAflow/`

Contains **ENEZAflow**, a reproducible Nextflow pipeline for 16S rRNA microbiome analysis. 

The pipeline processes raw sequencing reads through quality control, denoising, taxonomic assignment, and generation of analysis-ready outputs.

Detailed documentation is available in `ENEZAflow/README.md`.

### `Machine-Learning`

```text
Machine-Learning/
├── data/
├── figures/
├── model/
├── notebooks/
└── results/
```

This directory contains all resources related to machine learning model development, including:

* Processed datasets
* Model training scripts
* Trained models
* Performance evaluation results
* Figures

#### `notebooks/`

Contains seven well-documented Jupyter notebooks covering the complete machine learning workflow, from exploratory data analysis and preprocessing to model training, evaluation, and interpretation.

#### `model/`

Contains the trained and optimized machine learning models.

#### `figures/`

Contains figures generated during exploratory analysis, model evaluation, and interpretation.

#### `results/`

Contains model evaluation metrics, summaries, and exported outputs.

### `Results/`

Contains outputs generated from microbiome analyses, including diversity analyses, taxonomic summaries, and other downstream results.

### `Scripts/`

Contains utility scripts used throughout the project for data preprocessing, analysis, visualization, and workflow automation.
