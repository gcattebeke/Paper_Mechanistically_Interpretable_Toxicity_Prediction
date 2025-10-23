# High-Throughput Screening (HTS) Data

This directory contains the steps for extracting, filtering, and processing high-throughput screening data from the EPA's ToxCast/Tox21 programs via InvitroDB v4.2. Instructions and downloads to set up a local MySQL instance of InvitroDB, along with release notes, assay descriptions, summary files, and concentration–response plots, are available via EPA’s ToxCast data portal and the invitrodb v4.2 release space. 

The HTS data pipeline extracts bioactivity data from InvitroDB v4.2 and applies rigorous quality control filters to ensure high-quality datasets suitable for machine learning applications. The pipeline focuses on binary hit/no-hit classifications with high confidence labels. `input/TOX21_aeids_all.csv` contains the list of assay endpoint IDs (aeids) to be processed. The processing pipeline extracts data for these assays, applies multiple filtering steps, and outputs a final chemical-assay activity matrix in `output/all_assays_merged.csv`.

<br> 

## Directory Structure

```
HTS_data/
├── input/
│   └── TOX21_aeids_all.csv                             # list of assay endpoint IDs (aeid) to process
├── output/
│   └── all_assays_merged.csv                           # chemical-assay activity matrix after filtering
└── scripts/
    └── assay_database_extraction_and_selection.py      # processing pipeline
```

<br> 

## Pipeline Description

The processing pipeline `assay_database_extraction_and_selection.py` consists of five main steps:
1. **Database extraction**  
  Query InvitroDB v4.2 per aeid using an SQL template; apply QC and label logic (high‑confidence fit categories — Active: `41,42,37,38`; Inactive: `13,15`; require `hitc >= 0.90`); collapse replicates to a single per-chemical hit call via majority voting where appropriate.

2. **Pattern filtering**  
  Exclude non‑analytical assays matching configured patterns (defaults: `ch[12]`, `Followup`, `_RT_`, `AutoFluor`) to remove channel‑specific, follow‑up, time‑series, and autofluorescence assays.

3. **Viability filtering**  
  Identify chemicals flagged cytotoxic in paired viability assays (`hitc == 1`) and remove them from corresponding activity assays to avoid confounding by cell death.

4. **Class balance filtering**  
  Retain only assays with sufficient examples per class (default `MIN_CLASS_COUNT = 200` for both active and inactive); drop assays that fail this threshold.

5. **Matrix assembly**  
  Merge remaining assays into a wide matrix (rows = `dtxsid`, columns = assays), leaving missing values for untested combos and writing the final table to `output/all_assays_merged.csv`.
