
# ENEZAflow

*A Nextflow pipeline for reproducible 16S rRNA amplicon microbiome analysis.*

## Overview

**ENEZAflow** is a modular Nextflow pipeline developed to process Illumina 16S rRNA amplicon sequencing data from raw FASTQ files to high-quality amplicon sequence variants (ASVs) and downstream microbiome analyses. The pipeline emphasizes reproducibility, and scalability. The current pipeline has not been configured for use with containers but rather needs pre-installed `R`, `DADA2`, and `phyloseq` library.

The current implementation is optimized for paired-end Illumina sequencing data and uses the DADA2 algorithm for denoising and ASV inference.

---

## Features

* Quality assessment of raw sequencing reads
* Adapter and quality trimming
* DADA2-based denoising
* Paired-end read merging
* Chimera removal
* Taxonomic assignment using SILVA reference databases
* Construction of a Phyloseq object
* Export of ASV, taxonomy, and metadata tables
* MultiQC summary reports
* Reproducible execution using Nextflow

---

## Workflow

```text
Raw FASTQ
    |
    
Quality Control (FastQC)
    │
    
Read Trimming
    │
    
DADA2 Error Learning
    │
    
Dereplication
    │
    
ASV Inference
    │
    
Merge Paired Reads
    │
    
Chimera Removal
    │
    
Taxonomic Assignment (SILVA)
    │
    
Phyloseq Object
    |
    |- ASV Table
    |-- Taxonomy Table
    |-- Sample Metadata
    |-- Summary Reports
    |-- ML feature matrix
```

---

## Repository Structure

```text
ENEZAflow/
|-- bin/
|-- modules/
|--results/
|-- workflows/
|-- resources/
|-- main.nf
|-- nextflow.config
|-- README.md
```

---

## Requirements

* Nextflow ≥ 26
* Java 17 
* Base R v4.6.1
* Libraries (dada, optparse,phyloseq)
* Linux (recommended)

---

## Running the Pipeline
*To test the pipeline:
```bash
git clone

cd Eneza-Data-Science-Residential-training-2026/GUT_MICROBIOME_PRJCT_8/ENEZAflow


```
 & 
 Then run
```bash
nextflow run main.nf
```

- Make sure to edit the nextflow.config with appropriate file absolute paths for `input reads`, `metadata` and `silva ref sequences`.

---

## Input

The pipeline expects a path to paired-end FASTQ files and need to be in the form `*_1.fastq.gz` and `*_2.fastq.gz`
The pipeline will use the metadata to convert the input into a standard samplesheet.

Example:

| sample   | fastq_1                   | fastq_2                   |
| -------- | ------------------------- | ------------------------- |
| Sample01 | data/Sample01_1.fastq.gz | data/Sample01_2.fastq.gz |
| Sample02 | data/Sample02_1.fastq.gz | data/Sample02_2.fastq.gz |

---

## Output

The pipeline generates:

* Quality control reports
* Filtered reads
* ASV abundance table
* Taxonomic assignments
* Sequence table
* Phyloseq object
* MultiQC report
* Conversion of the abundance/asv tables into ML feature table
* Execution logs and reports

---

## Reference Databases

The SILVA reference databases are **not distributed with this repository** because of their size.

Download the required reference files separately and place them in the `resources/` directory before running the pipeline.

You can download the SILVA DB `_train_set.fa.gz` file from this [link](https://zenodo.org/records/4587955) and place it inside resources folder in the project dir.

---


---


## Authors

Developed as part of the **ENEZA Data Science Initiative Residential Training 2026**  for reproducible microbiome bioinformatics nextflow pipeline and downstream machine learning models for Type 2 Diabetes.
