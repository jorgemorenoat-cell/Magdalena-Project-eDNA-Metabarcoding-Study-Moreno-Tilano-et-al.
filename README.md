# Magdalena Project – eDNA Metabarcoding Study

## Overview

This repository contains all data, metadata, and scripts required to reproduce the analyses presented in Moreno-Tilano et al.(Preprint doi:10.22541/au.173262274.49229664/v1) The project focuses on environmental DNA (eDNA) metabarcoding to assess vertebrate diversity in the Magdalena River.

The repository includes processed outputs, taxonomic assignments, curated datasets, and scripts for bioinformatic processing and statistical analyses.

---

## Repository Structure

### 1. config/

This folder contains configuration files required for trimming, demultiplexing, and primer removal.

* `config_teleo.sh`: Parameter file used by the trimming pipeline (modifiable).
* `Corr_tags_SC22171_V05.csv`: Index file for Vert01 primer.
* `Corr_tags_SC22171_Teleo.csv`: Index file for Teleo primer.
* `all_samples.csv`: Sample metadata linking sample IDs to primers.

  * Samples starting with **SPY** correspond to this study.
  * Other samples belong to external projects but are included in the raw sequencing files.

---

### 2. outputs/

This folder contains processed data and intermediate results generated from the analyses.

Subfolders:

* `teleo/` and `vert01/`: Contain primer-specific outputs.

  * Taxonomic assignment tables.
  * `Taxa_counts.csv`: Taxon counts derived from subsampling analyses for sequencing depth evaluation.

Main files:

* `teleo_dada2RDP_...csv` and `Vert01_dada2RDP_...csv`: Raw taxonomic assignment tables from DADA2 processing.
* `Magdalena_data_2022_teleo_vert_dada2.csv`: Combined non-curated dataset after merging both markers and filtering by read counts.
* `List of vertebrates from all sources in the study area.csv`: Vertebrate records downloaded from GBIF for the study area.
* Excel file containing:

  * Curated dataset
  * Non-curated dataset
  * Simple verification dataset
  * Abbreviation tables

---

### 3. Scripts/

This folder contains all scripts required for bioinformatic processing and downstream analyses.

* `trimming.sh`: Main script for demultiplexing and primer trimming.
* `bash_functions.sh`: Auxiliary bash functions used by the trimming script.
* Script for subsampling and DADA2 workflow:

  * Generates subsamples of raw reads
  * Executes automated DADA2 pipeline
  * Produces taxonomic assignment and taxon counts
* `Dada2_automatico_bash_taxones`: Required for running the automated workflow.
* `Moreno-Tilano et al.R`: R script containing all statistical analyses and figure generation used in the manuscript.

---

### 4. Teleo Raw data/

### 5. Vert01 Raw data/

These folders contain raw sequencing data (compressed FASTQ files), including forward and reverse reads prior to demultiplexing.

---

### 6. utils/

This folder contains supporting datasets and metadata:

* `Faunal list.xlsx`: Compiled local faunal list used for taxonomic curation.
* Reference databases downloaded using CRABS.
* `Abreviaturasf.xlsx`:

  * Sheet 1: Sampling site coordinates
  * Sheet 2: Mapping between capsule codes and sample names used in the manuscript

---

## Reproducibility Workflow

To reproduce the analyses:

1. Start from raw data (`Teleo Raw data/`, `Vert01 Raw data/`)
2. Run trimming and demultiplexing:

   * `trimming.sh` (uses configuration files in `config/`)
3. Process sequences using DADA2:

   * The codes are present at the beginning of the main script.
4. Generate taxonomic assignments and ASV tables
5. Use curated datasets from `outputs/` for ecological analyses
6. Run statistical analyses and figure generation:

   * `Moreno-Tilano et al.R`

---

## Notes on Taxonomic Curation

Taxonomic curation was conducted manually following a structured decision framework described in the manuscript (Figure 2 and Figure S2). Curated and non-curated datasets are both provided to ensure transparency and reproducibility.

---

## Data Availability

All files included in this repository are sufficient to reproduce the results presented in the manuscript. Additional details on methods and analytical decisions are described in the associated publication.

---

## Contact

For questions regarding the dataset or scripts, please contact the corresponding author.
