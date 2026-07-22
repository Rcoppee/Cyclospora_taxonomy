# Cyclospora Population Genomics

**Genome-wide population genomics of *Cyclospora cayetanensis* using targeted hybridization capture**

---

## Overview

This repository contains the bioinformatic workflows and analysis scripts used in the manuscript:

> **Targeted genome capture reveals the global population structure of *Cyclospora cayetanensis* and enables molecular surveillance.**

The objective of this repository is to provide a fully reproducible workflow, from raw Illumina sequencing reads to genome-wide population genomic analyses and the identification of diagnostic molecular markers.

The repository includes scripts for:

- Whole-genome read mapping
- SNP discovery and variant filtering
- Population genomic analyses
- Identification of highly discriminatory genomic regions
- Generation of publication figures

---

## Repository structure

```text
Cyclospora_population_genomics/
│
├── assembling/
│   ├── Quality control
│   ├── Read trimming
│   ├── Read mapping
│   ├── BAM processing
│   └── Coverage analyses
│
├── variant_calling/
│   ├── SNP calling
│   ├── Variant filtering
│   ├── VCF merging
│   └── SNP matrix generation
│
├── population_structure/
│   ├── Phylogenetic analyses
│   ├── Principal Component Analysis (PCA)
│   ├── Discriminant Analysis of Principal Components (DAPC)
│   ├── Identity-by-Descent (IBD)
│   ├── Admixture analyses
│   └── Population assignment
│
├── marker_discovery/
│   ├── Genome-wide SNP screening
│   ├── Identification of highly informative loci
│   ├── Geographic marker discovery
│   └── Extraction of candidate regions for molecular assay design
│
├── figures/
│   Scripts used to reproduce all manuscript figures
│
└── README.md
```

---

# Workflow

The complete workflow consists of four major analytical steps.

## 1. Genome assembly and read mapping

Raw paired-end Illumina reads are processed and aligned against the *Cyclospora cayetanensis* reference genome (strain **NF1_C8**).

The workflow includes:

- FASTQ quality assessment
- Adapter trimming
- Read mapping
- BAM sorting and indexing
- Duplicate removal
- Mapping statistics
- Genome coverage calculation

### Output

- Sorted BAM files
- Mapping statistics
- Coverage metrics

---

## 2. Variant calling

Genome-wide SNPs are identified from aligned reads.

This workflow performs:

- SNP calling
- Quality filtering
- Removal of low-confidence variants
- VCF merging
- SNP matrix generation

### Output

- Filtered VCF files
- SNP matrices
- Genotype tables

---

## 3. Population genomic analyses

This directory contains all scripts used to reproduce the analyses presented in the manuscript.

Analyses include:

- Neighbor-Joining phylogenetic reconstruction
- Principal Component Analysis (PCA)
- Discriminant Analysis of Principal Components (DAPC)
- Population assignment analyses
- Identity-by-Descent (IBD) analyses
- Pairwise genetic distances
- Population differentiation analyses

These scripts also generate the publication-ready figures.

---

## 4. Marker discovery

Genome-wide SNP datasets were screened to identify genomic regions showing high discriminatory power between major geographic populations.

These analyses were used to:

- identify highly informative SNPs;
- extract candidate genomic regions;
- design molecular assays for:
  - Sanger sequencing;
  - High-Resolution Melting (HRM) qPCR.

---

# Input data

The workflows require:

- paired-end Illumina FASTQ files;
- *Cyclospora cayetanensis* reference genome (FASTA);
- reference annotation (optional).

Reference genome:

> *Cyclospora cayetanensis* strain **NF1_C8**

---

# Software requirements

The pipeline was developed using the following software:

| Software | Version |
|-----------|---------|
| BWA-MEM2 | XX |
| Samtools | XX |
| BCFtools | XX |
| GATK | XX |
| Mosdepth | XX |
| BEDTools | XX |
| FastQC | XX |
| fastp | XX |
| R | XX |
| Python | XX |

Required R packages are specified within each script.

---

# Reproducibility

Unless otherwise specified:

- default software parameters were used;
- identical filtering criteria were applied to all samples;
- random seeds were fixed whenever stochastic procedures were involved.

The scripts are organized to reproduce all analyses presented in the manuscript.

---

# Data availability

Raw sequencing reads generated during this study have been deposited in the **European Nucleotide Archive (ENA)**.

Accession numbers are provided in the manuscript.

---

# Citation

If you use this repository, please cite:

> Coppée R. *et al.*  
> **Targeted genome capture reveals the global population structure of *Cyclospora cayetanensis* and enables molecular surveillance.**

---

# Contact

**Romain Coppée**

French National Reference Centre for Malaria

French National Reference Centre for Cryptosporidiosis, Microsporidiosis and other Digestive Protozoa

Université de Rouen Normandie

Rouen, France

📧 romain.coppee@univ-rouen.fr *(replace if needed)*

---

# License

This project is distributed under the **MIT License**.
