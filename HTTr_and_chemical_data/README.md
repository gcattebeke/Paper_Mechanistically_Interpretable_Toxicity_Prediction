# High-Throughput Transcriptomics (HTTr) Data

This directory contains the complete pipeline for processing high-throughput transcriptomics data from three cell lines (MCF-7, U2OS, and HepRG), including differential expression analysis, pathway enrichment, and integration with chemical structures and HTS assay data.

<br> 

## Directory Structure

```
HTTr_data/
├── input/
│   ├── count_data/                    # Raw count data (see count_data/README.md)
│   └── metadata/                      # Sample metadata and chemical annotations
│       ├── File_S5_httr_mcf7_screen_qc_metrics.csv
│       ├── Sample_key.rds
│       ├── GSE284321_metadata.xlsx
│       ├── Chemical_dictionary.rds
│       ├── chemical_dtxsid_map.csv
│       ├── hWTv1_mcf7_probes.csv
│       ├── SMILES_data.csv
│       └── CCD-Batch-Search_2025-05-19_08_42_02.csv
├── output/
│   ├── aggregated_results_MCF7.feather
│   ├── aggregated_results_U2OS.feather
│   ├── aggregated_results_HepRG.feather
│   ├── signature_scores_all_databases.feather
│   └── chemical_httr_assay_aggregated.feather
└── scripts/
    ├── 01a_MCF7_DE_analysis.Rmd
    ├── 01b_U2OS_DE_analysis.Rmd
    ├── 01c_HepRG_DE_analysis.Rmd
    ├── 02_combining_and_activity_mapping.Rmd
    └── 03_chemical_intergration_and_data_aggregation.ipynb
```

## Pipeline Overview

The HTTr data processing pipeline consists of three main stages:

1. **Differential Expression Analysis (Scripts 01a-01c)**
Process raw count data for each cell line independently:
- Quality control filtering
- DESeq2 differential expression analysis
- Probe-to-gene aggregation
- Dose-response profiling

2. **Activity Mapping (Script 02)**
Combine results and compute biological signatures:
- Merge data from all three cell lines
- Calculate pathway activities (PROGENy)
- Calculate transcription factor activities (DoRothEA)
- Calculate gene set enrichment scores (MSigDB H + C2)

3. **Chemical Integration (Script 03)**
Integrate with chemical structures and HTS data:
- Generate chemical fingerprints (MACCS keys)
- Merge HTTr signatures with HTS assay data
- Aggregate dose-response curves into summary features
- Produce final integrated dataset