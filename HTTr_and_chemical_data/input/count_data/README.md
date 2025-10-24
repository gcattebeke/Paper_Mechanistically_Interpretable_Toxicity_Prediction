# High-Throughput Transcriptomics (HTTr) Count Data

This directory contains the raw count data files for high-throughput transcriptomics experiments across three cell lines.

## Data Sources

### MCF-7 (Human Breast Cancer Epithelial Cells)
- **Source**: NCBI Gene Expression Omnibus (GEO)
- **Accession**: [GSE272548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272548)
- **Platform**: L1000 gene expression profiling
- **Format**: Multiple CSV files with plate group-specific count data
- **File pattern**: `pg*_count_data.csv`

### U2OS (Human Bone Osteosarcoma Epithelial Cells)
- **Source**: Clowder Data Management System
- **URL**: [Clowder Dataset](https://clowder.edap-cluster.com/datasets/61147fefe4b0856fdc65639b#folderId=65959df8e4b063812d5c7706)
- **Platform**: L1000 gene expression profiling
- **Format**: RDS file with probe-level counts
- **File**: `Probe_counts.rds`

### HepRG (Human Hepatocyte Cell Line)
- **Source**: NCBI Gene Expression Omnibus (GEO)
- **Accession**: [GSE284321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284321)
- **Platform**: L1000 gene expression profiling
- **Format**: Multiple CSV files with plate group-specific count data
- **File pattern**: `pg*_count_data.csv`

## File Structure

```
count_data/
├── README.md                          # This file
├── pg1_count_data.csv                 # MCF-7 or HepRG plate group 1 counts
├── pg2_count_data.csv                 # MCF-7 or HepRG plate group 2 counts
├── ...                                # Additional plate groups
└── Probe_counts.rds                   # U2OS probe counts (R data format)
```

## Data Format

### CSV Files (MCF-7 and HepRG)
- **Rows**: Probe/gene identifiers
- **Columns**: Sample identifiers
- **Values**: Raw or normalized count values
- **First column**: Probe/gene names or IDs

### RDS File (U2OS)
- **Format**: R data frame serialized as RDS
- **Rows**: Probe identifiers
- **Columns**: Sample identifiers
- **Values**: Raw or normalized count values

## Download Instructions

### MCF-7 Data
1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272548
2. Navigate to "Download family" section
3. Download supplementary files containing count data
4. Extract CSV files to this directory

### U2OS Data
1. Visit: https://clowder.edap-cluster.com/datasets/61147fefe4b0856fdc65639b#folderId=65959df8e4b063812d5c7706
2. Log in or request access if required
3. Download `Probe_counts.rds` file
4. Place in this directory

### HepRG Data
1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284321
2. Navigate to "Download family" section
3. Download supplementary files containing count data
4. Extract CSV files to this directory

## Data Processing

These raw count files are processed by the differential expression analysis scripts:
- `01a_MCF7_DE_analysis.Rmd`: Processes MCF-7 count data
- `01b_U2OS_DE_analysis.Rmd`: Processes U2OS count data
- `01c_HepRG_DE_analysis.Rmd`: Processes HepRG count data

## File Size Notes

Count data files can be large (100s of MB to several GB). This directory is typically excluded from Git tracking via `.gitignore` and users should download the data directly from the source repositories.

## Quality Control

All count data undergoes quality control filtering during the differential expression analysis:
- Samples must pass QC flags provided in metadata
- Low-count genes/probes are filtered (mean count ≥ 5)
- Bad probes are excluded based on probe flags

## Related Files

- **Metadata**: Located in `../metadata/` directory
  - `File_S5_httr_mcf7_screen_qc_metrics.csv`: MCF-7 sample metadata and QC
  - `Sample_key.rds`: U2OS sample metadata and QC
  - `GSE284321_metadata.xlsx`: HepRG sample metadata
- **Probe annotations**: `../metadata/hWTv1_mcf7_probes.csv`

## Technical Specifications

### L1000 Platform
The L1000 platform measures the expression of approximately 978 "landmark" genes directly and infers the expression of an additional ~10,000 genes computationally. These experiments use the landmark genes for differential expression analysis.

### Count Normalization
- Count data may be pre-normalized or raw depending on the source
- DESeq2 normalization is applied during differential expression analysis
- Plate effects are controlled for in the statistical model

## Citation

If you use this data, please cite the original data sources:
- GSE272548: [Citation for MCF-7 study]
- GSE284321: [Citation for HepRG study]
- Clowder U2OS dataset: [Citation for U2OS study]

## Support

For issues with downloading or processing count data, please:
1. Check the original data repositories for file availability
2. Verify file paths in the R markdown analysis scripts
3. Open an issue on the main repository