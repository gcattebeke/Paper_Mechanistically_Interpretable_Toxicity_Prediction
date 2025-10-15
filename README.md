# Toxicogenomics Data Analysis Pipeline

This repository contains the complete data analysis pipeline and supplementary materials for the toxicogenomics research paper. The repository integrates high-throughput screening (HTS) data from ToxCast/Tox21 assays with high-throughput transcriptomics (HTTr) data to create a comprehensive toxicological dataset suitable for machine learning applications.

## Repository Structure

```
toxicogenomics_paper/
├── HTS_data/               # High-Throughput Screening data processing
│   ├── input/             # Raw assay endpoint identifiers
│   ├── output/            # Processed HTS activity data
│   └── scripts/           # Database extraction and filtering pipeline
│
└── HTTr_data/             # High-Throughput Transcriptomics data processing
    ├── input/             # Raw count data and metadata
    ├── output/            # Processed transcriptomic signatures
    └── scripts/           # Differential expression and integration pipeline
```

## Overview

This pipeline consists of two main components:

### 1. High-Throughput Screening (HTS) Data
- **Source**: InvitroDB v4.2 (ToxCast/Tox21 assays)
- **Purpose**: Extract and filter high-quality bioactivity data
- **Output**: Chemical-assay activity matrix with consensus hit calls

### 2. High-Throughput Transcriptomics (HTTr) Data
- **Sources**: 
  - MCF-7 cells: [GSE272548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272548)
  - U2OS cells: [Clowder Dataset](https://clowder.edap-cluster.com/datasets/61147fefe4b0856fdc65639b#folderId=65959df8e4b063812d5c7706)
  - HepRG cells: [GSE284321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284321)
- **Purpose**: Perform differential expression analysis and pathway enrichment
- **Output**: Integrated chemical-transcriptomic-assay dataset with aggregated features

## Requirements

### Python Dependencies
```
python >= 3.10.16
pandas >= 2.2.3
sqlalchemy >= 2.0.40
mysql-connector-python
rdkit
scipy
joblib
pyarrow
```

### R Dependencies
```
R >= 4.3.3
DESeq2 >= 1.42.1
dplyr >= 1.1.4
tidyverse >= 2.0.0
GSVA >= 1.50.5
decoupleR >= 2.8.0
dorothea >= 1.14.1
progeny >= 1.24.0
arrow >= 21.0.0.1
BiocParallel >= 1.36.0
```

### Database Requirements
- Access to InvitroDB v4.2 (MySQL database)
- Set environment variables:
  ```bash
  export DB_USER="your_username"
  export DB_PASS="your_password"
  export DB_HOST="localhost"
  export DB_PORT="3306"
  export DB_NAME="invitrodb_v4_2"
  ```

## Usage

### Step 1: Process HTS Data
```bash
cd HTS_data/scripts
python assay_database_extraction_and_selection.py
```
This script:
1. Extracts assay data from InvitroDB v4.2
2. Filters non-analytical file types
3. Removes cytotoxicity-confounded hits
4. Ensures class balance (≥200 samples per class)
5. Produces `all_assays_merged.csv` with chemical-assay activity matrix

### Step 2: Process HTTr Data

#### 2a. Differential Expression Analysis
Run the R markdown notebooks in order:
```r
# In RStudio or R console
rmarkdown::render("HTTr_data/scripts/01a_MCF7_DE_analysis.Rmd")
rmarkdown::render("HTTr_data/scripts/01b_U2OS_DE_analysis.Rmd")
rmarkdown::render("HTTr_data/scripts/01c_HepRG_DE_analysis.Rmd")
```

#### 2b. Combine and Map Activities
```r
rmarkdown::render("HTTr_data/scripts/02_combining_and_activity_mapping.Rmd")
```
This step:
- Merges results from all three cell lines
- Computes pathway activity scores (PROGENy)
- Computes transcription factor activities (DoRothEA)
- Computes gene set enrichment scores (MSigDB)

#### 2c. Final Integration
```bash
# Run the Jupyter notebook
jupyter notebook HTTr_data/scripts/03_chemical_intergration_and_data_aggregation.ipynb
```
This step:
- Generates chemical fingerprints (MACCS keys)
- Integrates HTTr signatures with HTS assay data
- Aggregates dose-response curves into summary features
- Produces final `chemical_httr_assay_aggregated.feather` file

## Output Files

### HTS Data
- `HTS_data/output/all_assays_merged.csv`: Chemical-by-assay binary activity matrix (dtxsid × assay_name)

### HTTr Data
- `HTTr_data/output/aggregated_results_MCF7.feather`: Gene-level log2 fold changes for MCF-7
- `HTTr_data/output/aggregated_results_U2OS.feather`: Gene-level log2 fold changes for U2OS
- `HTTr_data/output/aggregated_results_HepRG.feather`: Gene-level log2 fold changes for HepRG
- `HTTr_data/output/signature_scores_all_databases.feather`: Pathway and TF activity scores
- `HTTr_data/output/chemical_httr_assay_aggregated.feather`: Final integrated dataset

## Data Processing Methods

### HTS Pipeline
1. **Database Extraction**: Queries InvitroDB for high-confidence fit categories
2. **Pattern Filtering**: Removes channel-specific, follow-up, time-series, and autofluorescence assays
3. **Viability Filtering**: Excludes chemicals flagged as cytotoxic in paired viability assays
4. **Class Balance**: Retains only assays with ≥200 actives and ≥200 inactives
5. **Consensus Hit Calling**: Resolves replicates using majority voting (hitc ≥ 0.90)

### HTTr Pipeline
1. **Quality Control**: Filters samples based on provided QC flags
2. **Differential Expression**: DESeq2 with plate and dose as covariates
3. **Gene-Level Aggregation**: Maximum absolute fold change across probes
4. **Pathway Analysis**: 
   - PROGENy: 14 cancer-relevant pathways
   - DoRothEA: Transcription factor activities (confidence A/B/C)
   - MSigDB: Hallmark + C2 canonical pathways (ssGSEA)
5. **Dose-Response Aggregation**: 
   - Maximum robust value (median of top-2)
   - Dose at maximum response
   - Area under curve for negative responses

## Citation

If you use this code or data in your research, please cite:

[Paper citation to be added upon publication]

## License

[License information to be added]

## Contact

For questions or issues, please open an issue on this repository or contact the corresponding author.

## Acknowledgments

This work uses data from:
- EPA's ToxCast program and the Tox21 consortium
- Gene Expression Omnibus (GEO) datasets: GSE272548, GSE284321
- The Clowder data management platform
