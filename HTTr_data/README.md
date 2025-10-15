# High-Throughput Transcriptomics (HTTr) Data

This directory contains the complete pipeline for processing high-throughput transcriptomics data from three cell lines (MCF-7, U2OS, and HepRG), including differential expression analysis, pathway enrichment, and integration with chemical structures and HTS assay data.

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

### Stage 1: Differential Expression Analysis (Scripts 01a-01c)
Process raw count data for each cell line independently:
- Quality control filtering
- DESeq2 differential expression analysis
- Probe-to-gene aggregation
- Dose-response profiling

### Stage 2: Activity Mapping (Script 02)
Combine results and compute biological signatures:
- Merge data from all three cell lines
- Calculate pathway activities (PROGENy)
- Calculate transcription factor activities (DoRothEA)
- Calculate gene set enrichment scores (MSigDB H + C2)

### Stage 3: Chemical Integration (Script 03)
Integrate with chemical structures and HTS data:
- Generate chemical fingerprints (MACCS keys)
- Merge HTTr signatures with HTS assay data
- Aggregate dose-response curves into summary features
- Produce final integrated dataset

## Data Sources

### MCF-7 (Human Breast Cancer Epithelial Cells)
Count files were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272548

### U2OS (Human Bone Osteosarcoma Epithelial Cells)
Count files were downloaded from https://clowder.edap-cluster.com/datasets/61147fefe4b0856fdc65639b#folderId=65959df8e4b063812d5c7706

### HepRG (Human Hepatocyte Cell Line)
Count files were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284321

## Input Files

### Count Data
Located in `input/count_data/`. See `input/count_data/README.md` for details on downloading and formatting.

### Metadata Files

#### MCF-7 Metadata
- **File**: `File_S5_httr_mcf7_screen_qc_metrics.csv`
- **Contents**: Sample IDs, dose levels, plate groups, QC flags
- **Key columns**: `sample_id`, `epa_sample_id`, `dose_level`, `pg_id`, `qc_flag`, `dtxsid`

#### U2OS Metadata
- **File**: `Sample_key.rds` (R data format)
- **Contents**: Sample IDs, chemical IDs, dose levels, QC flags
- **Key columns**: `sample_id`, `chem_id`, `dose_level`, `qc_flag`

#### HepRG Metadata
- **File**: `GSE284321_metadata.xlsx`
- **Contents**: Sample information from GEO submission
- **Key columns**: `chemical sample ID`, `chemical name`, `dose level`

#### Chemical Annotations
- **Chemical_dictionary.rds**: Maps chemical IDs to DTXSIDs for U2OS
- **chemical_dtxsid_map.csv**: Manual curation of chemical names to DTXSIDs for HepRG
- **CCD-Batch-Search_2025-05-19_08_42_02.csv**: Batch search results from EPA CompTox Dashboard
- **SMILES_data.csv**: Chemical structures in SMILES format for fingerprint generation

#### Probe Annotations
- **hWTv1_mcf7_probes.csv**: L1000 probe annotations
- **Columns**: `Probe_Name`, `Gene_Symbol`, `Entrez_ID`, `Probe_Flag`

## Output Files

### Cell Line-Specific Results (Stage 1)

#### `aggregated_results_MCF7.feather`
Gene-level differential expression results for MCF-7 cell line.

**Format**: Apache Feather (binary columnar format)

**Key columns**:
- `sample_dose`: Sample ID and dose level combination
- `lfc_[GENE]`: Log2 fold change for each gene
- `padj_[GENE]`: Adjusted p-value for each gene

**Dimensions**: ~2,000 sample-dose combinations × ~1,000 genes

#### `aggregated_results_U2OS.feather`
Gene-level differential expression results for U2OS cell line (same format as MCF-7).

#### `aggregated_results_HepRG.feather`
Gene-level differential expression results for HepRG cell line (same format as MCF-7).

### Combined Signatures (Stage 2)

#### `signature_scores_all_databases.feather`
Pathway and transcription factor activity scores across all cell lines.

**Format**: Apache Feather

**Key columns**:
- `epa_sample_id`: Sample identifier
- `dose_level`: Dose level tested
- `cell_type`: MCF7, U2OS, or HepRG
- `outcome_id`: DTXSID chemical identifier
- `sample_dose_cell`: Unique identifier (sample_dose_cell format)
- **MSigDB signatures** (5,000+ columns): ssGSEA normalized enrichment scores
- **DoRothEA TFs** (100+ columns): Transcription factor activity scores
- **PROGENy pathways** (14 columns): Cancer pathway activity z-scores

**Dimensions**: ~5,000 treatments × ~5,200 features

### Final Integrated Dataset (Stage 3)

#### `chemical_httr_assay_aggregated.feather`
Complete integrated dataset ready for machine learning applications.

**Format**: Apache Feather

**Key columns**:
- `outcome_id`: DTXSID chemical identifier
- `epa_sample_id`: Sample identifier
- `cell_type`: Cell line
- `chnm`: Chemical name
- **Chemical features** (~167 columns):
  - `MACCS_1` to `MACCS_167`: MACCS structural fingerprints
- **HTTr aggregate features** (~15,600 columns):
  - `[FEATURE]_max`: Robust maximum value (median of top-2 absolute values)
  - `[FEATURE]_dose_at_max`: Normalized dose at maximum response
  - `[FEATURE]_AUC_neg`: Area under curve for negative responses
- **HTS assay features** (~100 columns):
  - Binary activity calls for ToxCast/Tox21 assays

**Dimensions**: ~1,000 unique chemical-cell combinations × ~16,000 features

## Script Descriptions

### `01a_MCF7_DE_analysis.Rmd`
**Purpose**: Differential expression analysis for MCF-7 cell line

**Key steps**:
1. Load count data and metadata
2. Filter samples (QC pass, exclude controls and reference chemicals)
3. Iterate through plate groups
4. For each test chemical vs. vehicle control:
   - Construct DESeq2 model with plate and dose as covariates
   - Filter low-count genes (mean ≥ 5)
   - Perform differential expression with lfcShrink
   - Aggregate probe-level to gene-level (max absolute fold change)
5. Output wide-format matrix with sample-dose as rows

**Runtime**: ~2-4 hours (with 40 cores)

### `01b_U2OS_DE_analysis.Rmd`
**Purpose**: Differential expression analysis for U2OS cell line

**Differences from 01a**:
- Uses `Sample_key.rds` for metadata
- Single count data file (`Probe_counts.rds`) instead of plate-grouped files
- Otherwise identical methodology

**Runtime**: ~2-4 hours (with 40 cores)

### `01c_HepRG_DE_analysis.Rmd`
**Purpose**: Differential expression analysis for HepRG cell line

**Differences from 01a**:
- Uses Excel metadata file with complex chemical name mappings
- Requires manual curation and batch search for DTXSID mapping
- Plate groups extracted from file names

**Runtime**: ~2-4 hours (with 40 cores)

### `02_combining_and_activity_mapping.Rmd`
**Purpose**: Combine cell lines and compute biological activity signatures

**Key steps**:
1. **Load and merge DE results** from all three cell lines
2. **Annotate with DTXSIDs** using chemical dictionaries
3. **Prepare expression matrix** for activity mapping:
   - Transpose to gene × sample format
   - Filter genes with <95% finite values
   - Replace NAs with zeros
4. **Compute MSigDB enrichment** (ssGSEA):
   - Hallmark gene sets
   - C2 canonical pathways
   - Handle UP/DN paired signatures
5. **Compute TF activities** (DoRothEA + decoupleR):
   - Confidence A/B/C transcription factors
   - Consensus of ULM, MLM, and WSUM methods
6. **Compute pathway activities** (PROGENy):
   - 14 cancer-relevant pathways
   - Z-scored normalized activities
7. **Join all signatures** into single data frame

**Runtime**: ~30-60 minutes (with 40 cores)

### `03_chemical_intergration_and_data_aggregation.ipynb`
**Purpose**: Generate chemical fingerprints and integrate with HTS data

**Key steps**:
1. **Generate MACCS fingerprints**:
   - Parse SMILES strings using RDKit
   - Compute 167-bit MACCS keys
   - Handle parsing failures gracefully
2. **Merge datasets**:
   - Join chemical features with HTTr signatures
   - Join with HTS assay activity matrix
   - Remove columns with no variance
3. **Aggregate dose-response curves**:
   - Group by chemical-sample-cell combinations
   - For each feature, compute:
     - **Robust maximum**: Median of top-2 absolute values (reduces outlier influence)
     - **Dose at maximum**: Normalized dose (0-1) where maximum occurs
     - **Negative AUC**: Area under curve for negative responses only
4. **Output final dataset** as Feather file

**Runtime**: ~15-30 minutes (with 60 cores)

## Usage

### Prerequisites

#### R Environment
```r
install.packages(c("readr", "dplyr", "tidyverse", "gtools", "arrow", "stringr"))
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "BiocParallel", "GSVA"))
BiocManager::install(c("decoupleR", "dorothea", "progeny", "msigdbr"))
```

#### Python Environment
```bash
pip install pandas numpy scipy joblib pyarrow rdkit jupyter
```

### Running the Pipeline

#### Step 1: Differential Expression (Sequential)
```r
# In RStudio or R console
setwd("HTTr_data/scripts")

# MCF-7 analysis
rmarkdown::render("01a_MCF7_DE_analysis.Rmd")

# U2OS analysis
rmarkdown::render("01b_U2OS_DE_analysis.Rmd")

# HepRG analysis
rmarkdown::render("01c_HepRG_DE_analysis.Rmd")
```

**Note**: These must be run in order but are independent of each other.

#### Step 2: Activity Mapping
```r
rmarkdown::render("02_combining_and_activity_mapping.Rmd")
```

**Note**: Requires all three Stage 1 outputs to be complete.

#### Step 3: Chemical Integration
```bash
cd HTTr_data/scripts
jupyter notebook 03_chemical_intergration_and_data_aggregation.ipynb
# Run all cells
```

**Note**: Requires HTS data (`../../HTS_data/output/all_assays_merged.csv`) and Stage 2 output.

## Computational Requirements

### Memory
- **Stage 1**: 32-64 GB RAM (depending on cell line size)
- **Stage 2**: 64-128 GB RAM (large expression matrices)
- **Stage 3**: 32-64 GB RAM

### CPU
- Scripts are parallelized using `BiocParallel` (R) and `joblib` (Python)
- Recommended: 20-60 cores for optimal performance
- Can run on fewer cores by adjusting `workers` parameter:
  ```r
  register(MulticoreParam(workers = 10))  # R
  ```
  ```python
  Parallel(n_jobs=10)  # Python
  ```

### Storage
- **Input data**: 5-10 GB (count files)
- **Output data**: 2-5 GB (compressed Feather format)
- **Temporary files**: Minimal (R manages in memory)

## Methodological Details

### Differential Expression Model
```r
design = ~ plate_id + dose_level
```
- **plate_id**: Controls for plate-to-plate technical variation
- **dose_level**: Factor with levels {0, 1, 2, 3, ...} where 0 = vehicle control
- **Contrasts**: Each dose level vs. dose level 0 (vehicle)

### Gene Aggregation Strategy
From paper: *"DESeq2 probe-level L2FC results for each test chemical were summarized as a matrix in which the rows and columns correspond to the treatment concentrations and probes, respectively. For each probe-level L2FC matrix, L2FC values were aggregated to the gene level by using the highest magnitude fold change for any associated probe in either direction."*

**Implementation**:
```r
group_by(Gene_Symbol) %>%
slice_max(order_by = abs(log2FoldChange), with_ties = FALSE)
```

### Dose-Response Aggregation

#### Robust Maximum
Instead of simple maximum, uses median of top-2 absolute values to reduce outlier sensitivity:
```python
k = min(2, y.size)
top_vals = y[np.argsort(-np.abs(y))[:k]]
max_robust = np.median(top_vals)
```

#### Dose at Maximum
Normalized dose (0-1 scale) where absolute maximum occurs:
```python
dose_norm = dose / dose.max()
dose_at = dose_norm[np.argmax(np.abs(vals))]
```

#### Negative AUC
Captures sustained negative responses (gene repression, pathway inhibition):
```python
y_neg = np.where(y < 0, y, 0.0)
auc_neg = np.trapz(y_neg, x=dose_norm)
```

## Quality Control

### Sample-Level QC
- QC flags must be "OK" (provided in metadata)
- Exclude control samples: "QC sample", "reference chemical", "untreated control"
- Include only: "test sample", "vehicle control"

### Gene-Level QC
- Filter genes with mean count < 5 (DESeq2 recommendation)
- Exclude probes with "bad" probe flags
- Require ≥95% finite values across samples for activity mapping

### Chemical-Level QC
- SMILES parsing must succeed for MACCS fingerprint generation
- Failed SMILES are logged and excluded from final dataset
- Manual curation of chemical names to DTXSIDs for HepRG

## Data Formats

### Feather Format
All output files use Apache Feather format:
- **Advantages**: Fast read/write, compressed, cross-platform (R/Python)
- **Reading in R**: `arrow::read_feather()`
- **Reading in Python**: `pd.read_feather()`

### Wide vs. Long Format
- **Stage 1 outputs**: Wide format (sample-dose × genes)
- **Stage 2 output**: Wide format (sample-dose-cell × signatures)
- **Stage 3 output**: Wide format (chemical-cell × all features)

Wide format is preferred for machine learning applications.

## Troubleshooting

### Memory Issues
```r
Error: cannot allocate vector of size X GB
```
**Solutions**:
1. Reduce number of workers: `MulticoreParam(workers = 10)`
2. Process cell lines on separate machines
3. Increase system swap space
4. Use high-memory compute node

### Missing DTXSIDs
```r
Warning: X samples have missing outcome_id
```
**Solutions**:
1. Check chemical dictionaries are complete
2. Run EPA CompTox batch search for missing names
3. Update mapping files in metadata directory

### Failed SMILES Parsing
```python
Warning: could not parse SMILES 'XXXXX'
```
**Solutions**:
1. Check SMILES syntax in `SMILES_data.csv`
2. Update to correct SMILES from PubChem or CompTox
3. Remove problematic chemicals if not critical

### Slow Performance
**Solutions**:
1. Increase parallel workers
2. Use SSD for temporary files
3. Ensure adequate RAM to avoid swapping
4. Check network speed if database is remote

## References

### Data Sources
- GEO GSE272548: MCF-7 HTTr screen
- GEO GSE284321: HepRG HTTr screen
- Clowder: U2OS HTTr screen
- EPA CompTox Dashboard: Chemical annotations

### Methods
- **DESeq2**: Love, M.I., et al. (2014) Genome Biology
- **PROGENy**: Schubert, M., et al. (2018) Nature Communications
- **DoRothEA**: Garcia-Alonso, L., et al. (2019) Genome Research
- **GSVA**: Hänzelmann, S., et al. (2013) BMC Bioinformatics
- **decoupleR**: Badia-i-Mompel, P., et al. (2022) Bioinformatics Advances

## Contact

For questions about the HTTr data pipeline, please open an issue on the main repository.