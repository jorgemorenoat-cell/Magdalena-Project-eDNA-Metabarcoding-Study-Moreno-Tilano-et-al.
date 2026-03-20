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
### 4. Raw sequencing data (Zenodo)

Raw sequencing data are not stored in this repository due to file size constraints and are available in Zenodo:

DOI: 10.22541/au.173262274.49229664/v1

Two types of raw data are provided:

#### 4.1. Original Illumina output (non-demultiplexed)

These files correspond to the direct output from the sequencing platform (paired-end FASTQ files) prior to demultiplexing and primer trimming.

The dataset includes four sequencing runs identified by the suffix “AX” in the file names:

- **A1 and A5** → Vert01 primer (general vertebrates)
- **A2 and A7** → Teleo primer (fish-specific)

File naming follows the Illumina convention:

- `*_R1.fastq.gz`: forward reads  
- `*_R2.fastq.gz`: reverse reads  

Example:
- `VL324___MB0423A2___R1.fastq.gz`
- `VL324___MB0423A2___R2.fastq.gz`

These files require demultiplexing and primer trimming before downstream analyses.

#### 4.2. Demultiplexed and primer-trimmed data

Additionally, we provide preprocessed raw data organized into:

- `Teleo Raw data/`
- `Vert01 Raw data/`

These compressed folders contain reads already:

- demultiplexed by sample  
- trimmed for primer sequences  

These files can be used directly as input for DADA2 and downstream analyses without running the trimming pipeline.

#### Relationship with this repository

- Original Illumina data → require processing using `Scripts/trimming.sh`
- Demultiplexed data → can be used directly for downstream analyses
- Sample identities and indices → `config/` folder
- Processed outputs → `outputs/` folder

Users must download the raw data from Zenodo and place them in the expected directory structure before running the pipeline.
---

### 5. utils/

This folder contains supporting datasets and metadata:

* `Faunal list.xlsx`: Compiled local faunal list used for taxonomic curation.
* Reference databases downloaded using CRABS.
* `Abreviaturasf.xlsx`:

  * Sheet 1: Sampling site coordinates
  * Sheet 2: Mapping between capsule codes and sample names used in the manuscript

---
## Reproducibility Workflow

Two alternative workflows are available depending on the starting point of the analysis:

### Option 1: Full pipeline (from original raw data)

1. Download original Illumina raw data from Zenodo  
2. Run demultiplexing and primer trimming:

   * `trimming.sh` (uses configuration files in `config/`)

3. Process sequences using DADA2:

   * Code available at the beginning of `Moreno-Tilano et al.R`

4. Generate ASV tables and taxonomic assignments  
5. Use curated datasets from `outputs/` for ecological analyses  
6. Run statistical analyses and figure generation:

   * `Moreno-Tilano et al.R`

---

### Option 2: Reduced pipeline (from demultiplexed data)

1. Download preprocessed data (`Teleo Raw data/`, `Vert01 Raw data/`) from Zenodo  
2. Skip demultiplexing and primer trimming  
3. Process sequences directly using DADA2:

   * Code available at the beginning of `Moreno-Tilano et al.R`

4. Generate ASV tables and taxonomic assignments  
5. Use curated datasets from `outputs/` for ecological analyses  
6. Run statistical analyses and figure generation  

---

### Important note

The full pipeline (Option 1) represents the complete and fully reproducible workflow starting from raw sequencing output. The reduced pipeline (Option 2) is provided for convenience and faster reproducibility of downstream analyses.

---

## Notes on Taxonomic Curation

Taxonomic curation was conducted manually following a structured decision framework described in the manuscript (Figure 2 and Figure S2). Curated and non-curated datasets are both provided to ensure transparency and reproducibility.

---

## Data Availability

All files included in this repository are sufficient to reproduce the results presented in the manuscript. Additional details on methods and analytical decisions are described in the associated publication.

---

## Contact

For questions regarding the dataset or scripts, please contact the corresponding author.
