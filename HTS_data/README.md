# High-Throughput Screening (HTS) Data

This directory contains the steps for extracting, filtering, and processing high-throughput screening data from the EPA's ToxCast/Tox21 programs via InvitroDB v4.2. Instructions and downloads to set up a local MySQL instance of InvitroDB, along with release notes, assay descriptions, summary files, and concentration–response plots, are available via EPA’s ToxCast data portal and the invitrodb v4.2 release space. 

The HTS data pipeline extracts bioactivity data from InvitroDB v4.2 and applies rigorous quality control filters to ensure high-quality datasets suitable for machine learning applications. The pipeline focuses on binary hit/no-hit classifications with high confidence labels.

`input/TOX21_aeids_all.csv` contains the list of assay endpoint IDs (aeids) to be processed. The processing pipeline extracts data for these assays, applies multiple filtering steps, and outputs a final chemical-assay activity matrix in `output/all_assays_merged.csv`.

<br> 

### Directory Structure

```
HTS_data/
├── input/
│   └── TOX21_aeids_all.csv          # list of assay endpoint IDs (aeid) to process
├── output/
│   └── all_assays_merged.csv        # chemical-assay activity matrix after filtering
└── scripts/
    └── assay_database_extraction_and_selection.py  # processing pipeline
```

<br> 

## Pipeline Description

The processing pipeline (`assay_database_extraction_and_selection.py`) consists of five main steps:

### 1. Database Extraction
Queries InvitroDB v4.2 for each assay endpoint using SQL template with:
- **High-confidence fit categories**:
  - Active: fitc ∈ {41, 42, 37, 38}
  - Inactive: fitc ∈ {13, 15}
- **Hit call threshold**: hitc ≥ 0.90 for active classification
- **Consensus logic**: Resolves replicate measurements via majority voting

### 2. Pattern Filtering
Removes non-analytical assay types:
- **Channel-specific data** (`ch1`, `ch2`): Only ratio files provide normalized data
- **Follow-up assays** (`Followup`): Limited chemical coverage
- **Time-series data** (`_RT_`): Incompatible with endpoint analysis
- **Autofluorescence** (`AutoFluor`): Technical artifacts, not biological activity

### 3. Viability Filtering
Excludes compounds showing cytotoxicity in paired viability assays:
- Identifies chemicals flagged as toxic (hitc = 1) in viability assays
- Removes these chemicals from corresponding activity assays
- Prevents confounding of biological activity with cell death

### 4. Class Balance Filtering
Ensures sufficient samples for machine learning:
- Minimum 200 active samples (hitc = 1)
- Minimum 200 inactive samples (hitc = 0)
- Assays not meeting criteria are excluded

### 5. Matrix Assembly
Combines all filtered assays into a single wide-format matrix:
- Rows indexed by dtxsid (DSSTox Substance Identifier)
- Columns represent individual assays
- Missing values indicate untested chemical-assay combinations

## Usage

### Prerequisites
1. Access to InvitroDB v4.2 MySQL database
2. Set environment variables:
   ```bash
   export DB_USER="your_username"
   export DB_PASS="your_password"
   export DB_HOST="localhost"  # or remote host
   export DB_PORT="3306"
   export DB_NAME="invitrodb_v4_2"
   ```

### Installation
```bash
pip install pandas>=2.2.3 sqlalchemy>=2.0.40 mysql-connector-python
```

### Running the Pipeline
```bash
cd HTS_data/scripts
python assay_database_extraction_and_selection.py
```

### Expected Runtime
- **Database extraction**: 5-10 minutes (depending on network speed)
- **Filtering steps**: 2-5 minutes
- **Total runtime**: ~10-15 minutes

### Output
The script will print progress information:
```
Database extraction: 150 extracted, 0 skipped
Pattern filtering: 20/150 files excluded
Viability filtering: 45 assays processed
Class balance filtering: 25 assays excluded
Final dataset: 105 assays retained
```

## Configuration

### Filtering Thresholds
Edit these constants in `assay_database_extraction_and_selection.py`:

```python
# Fit categories (high-confidence labels)
FITC_STRONG_ACTIVES = (41, 42, 37, 38)
FITC_STRONG_INACTIVES = (13, 15)

# Hit call threshold for activity
HITC_ACTIVE_THRESHOLD = 0.90

# Minimum samples per class
MIN_CLASS_COUNT = 200

# Patterns to exclude
EXCLUDE_PATTERNS = [
    r'ch[12]',      # Channel-specific
    r'Followup',    # Follow-up assays
    r'_RT_',        # Time-series
    r'AutoFluor',   # Autofluorescence
]
```

## Temporary Files

The pipeline creates temporary processing directories (automatically cleaned up):
```
.temp_processing/
├── 00_raw_extracts/           # Initial database extracts
├── 01_pattern_filtered/       # After pattern filtering
├── 02_viability_filtered/     # After viability filtering
└── 03_class_balanced/         # After class balance filtering
```

These are removed after successful completion.

## Data Quality

### Reliability Features
- **Consensus hit calling**: Aggregates replicates to reduce noise
- **High-confidence labels**: Uses strict fit category criteria
- **Viability filtering**: Removes cytotoxicity confounders
- **Class balance**: Ensures sufficient data for both classes

### Quality Control Flags
The pipeline automatically excludes:
- Samples with QC flags (mc6_mthd_id ≠ 20)
- Low-confidence fit categories
- Assays with imbalanced or insufficient data

## Technical Notes

### Database Schema
The SQL query joins multiple InvitroDB tables:
- `mc5`: Level 5 modeling results (hit calls)
- `mc4`: Level 4 modeling results (curve fitting)
- `mc6`: Level 6 flags (quality control)
- `sample`: Sample identifiers
- `chemical`: Chemical annotations

### Performance
- Uses SQLAlchemy for efficient database connections
- Parallel processing can be added using `joblib` or `multiprocessing`
- Connection pooling via `pool_pre_ping=True`

## Troubleshooting

### Common Issues

**Database Connection Error**
```
sqlalchemy.exc.OperationalError: (2003, "Can't connect to MySQL server")
```
**Solution**: Check that DB_HOST, DB_PORT, and credentials are correct

**Empty Output**
```
Class balance filtering: 150 assays excluded
Final dataset: 0 assays retained
```
**Solution**: Lower `MIN_CLASS_COUNT` threshold or check input aeid list

**Missing Viability Assays**
```
Viability filtering: 0 assays processed
```
**Solution**: Normal if no viability assays in input list; filtering is optional

## References

- InvitroDB v4.2: [EPA CompTox Dashboard](https://www.epa.gov/chemical-research/exploring-toxcast-data)
- ToxCast/Tox21: [EPA ToxCast Program](https://www.epa.gov/chemical-research/toxicity-forecasting)
- Data Dictionary: [ToxCast Data Dictionary](https://www.epa.gov/chemical-research/toxcast-data-dictionary)

## Contact

For questions about the HTS data pipeline, please open an issue on the main repository.
