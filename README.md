# Mechanistically Interpretable Toxicity Prediction Through Multimodal Integration of Structure and Transcriptomics

![Workflow Overview](extra/overview_machine_learning_framework.png)

## Overview

This repository contains the complete computational pipeline for predicting chemical toxicity endpoints by integrating molecular structure fingerprints and dose-dependent transcriptomic signatures. The approach combines high-throughput transcriptomics (HTTr) data from three human cell lines (MCF7, U2OS, HepRG) with chemical structure information (MACCS fingerprints) to train interpretable machine learning models for 41 curated Tox21 assay endpoints.

**Key Innovation**: This framework leverages transcriptomic dose-response profiles (gene expression changes across concentration gradients) as mechanistic features, aggregated through pathway analysis (PROGENy, MSigDB, DoRothEA) to capture biological mechanisms underlying chemical toxicity.

---

## Repository Structure

```
toxicogenomics_paper/
├── analyses/               # Sequential analysis pipeline scripts
│   ├── 01_assay_database_extraction_and_selection.py
│   ├── 02a_MCF7_DE_analysis.Rmd
│   ├── 02b_U2OS_DE_analysis.Rmd
│   ├── 02c_HepRG_DE_analysis.Rmd
│   ├── 03_combining_and_activity_mapping.Rmd
│   ├── 04_chemical_intergration_and_data_aggregation.ipynb
│   ├── 05_modelling_and_evaluation.py
│   └── 06_results_analyses.ipynb
├── data/                   # Input data and intermediate outputs
│   ├── HTS/               # High-throughput screening assay data
│   ├── HTTr/              # Transcriptomic data and metadata
│   └── chemical/          # Chemical structure and fingerprint data
├── output/                # Model results and performance metrics
│   └── [41 assay directories with model outputs]
└── extra/                 # Supporting figures and documentation
```

---

## Pipeline Overview

### 1. Assay Database Extraction (`01_assay_database_extraction_and_selection.py`)

**Purpose**: Extract and curate high-throughput screening (HTS) data from InvitroDB v4.2 for Tox21 assays.

**Key Steps**:
- Queries InvitroDB v4.2 for 41 Tox21 assay endpoints (ACE inhibition, nuclear receptor activity, etc.)
- Applies consensus hit-calling logic across technical replicates
- Filters assays based on:
  - Analytical quality (removes non-endpoint data like channel-specific or time-series files)
  - Viability interference (excludes chemicals flagged as cytotoxic in paired viability assays)
  - Class balance (minimum 200 samples per class for robust modeling)
- Outputs: `all_assays_merged.csv` (chemical × assay matrix)

**Technologies**: Python 3.10, SQLAlchemy, pandas

---

### 2. Differential Expression Analysis (`02a-c_*_DE_analysis.Rmd`)

**Purpose**: Perform dose-dependent differential gene expression analysis for three cell lines using HTTr data.

**Data Sources**:
- **MCF7**: Breast cancer cell line (GSE272548)
- **U2OS**: Osteosarcoma cell line (Clowder repository)
- **HepRG**: Hepatocyte-like cell line (GSE284321)

**Key Steps**:
1. Load count data and metadata for each cell line
2. Run DESeq2 analysis with design: `~ plate_id + dose_level`
3. Filter genes with mean counts ≥ 5 across samples
4. Perform pairwise comparisons (each dose vs. vehicle control)
5. Apply lfcShrink with `ashr` method for regularized log2 fold changes
6. Aggregate probe-level results to gene-level (highest magnitude fold change)
7. Output: Wide-format matrices with samples × genes (log2FC and adjusted p-values)

**Technologies**: R 4.3.3, DESeq2, BiocParallel (40 cores)

---

### 3. Transcriptomic Activity Mapping (`03_combining_and_activity_mapping.Rmd`)

**Purpose**: Transform gene-level fold changes into interpretable pathway and transcription factor activity scores.

**Feature Engineering Methods**:
1. **MSigDB Hallmark + C2 Collections** (ssGSEA):
   - Single-sample gene set enrichment for canonical pathways
   - UP/DN signature pairs combined (UP - DN) for directionality
   
2. **DoRothEA Transcription Factor Activity** (decoupleR):
   - Consensus estimation across ULM, MLM, and weighted sum methods
   - High-confidence TF-target interactions (confidence A/B/C)
   
3. **PROGENy Pathway Activity**:
   - Quantifies activity in 14 cancer-related signaling pathways
   - Z-score normalized for cross-sample comparability

**Output**: `signature_scores_all_databases.feather` (samples × ~6,000 pathway/TF features)

**Technologies**: R 4.3.3, GSVA, decoupleR, progeny, msigdbr

---

### 4. Chemical Integration & Data Aggregation (`04_chemical_intergration_and_data_aggregation.ipynb`)

**Purpose**: Integrate chemical structure fingerprints with dose-aggregated transcriptomic features.

**Key Steps**:
1. **Chemical Structure Processing**:
   - Generate MACCS 166-bit molecular fingerprints from SMILES
   - Capture structural features (rings, functional groups, bond patterns)

2. **Dose-Response Aggregation** (per chemical-cell line pair):
   - **Maximum Response**: Robust max (median of top-2 absolute values)
   - **Dose at Maximum**: Normalized dose level where max occurs
   - **AUC Negative**: Area under curve for negative (repressive) responses
   
3. **Data Merging**:
   - Combine chemical fingerprints + aggregated transcriptomics + HTS labels
   - Remove constant features and samples with missing fingerprints

**Output**: `chemical_httr_assay_aggregated.feather` (chemicals × [166 MACCS + ~18,000 HTTr features + 41 assay labels])

**Technologies**: Python, RDKit, scikit-learn, pandas

---

### 5. Model Training & Evaluation (`05_modelling_and_evaluation.py`)

**Purpose**: Train XGBoost classifiers with nested cross-validation and SHAP interpretability for all 41 assays.

**Machine Learning Pipeline**:

**Cross-Validation Strategy**:
- 3-fold Stratified Group K-Fold (ensures same chemical never in train/test simultaneously)
- Nested CV with inner 3-fold for hyperparameter tuning

**Feature Selection**:
- **Boruta algorithm** on HTTr features (stability-based across 5 random undersamples)
- Retains features selected in ≥60% of runs for robustness
- Chemical fingerprints (MACCS) always included

**Class Imbalance Handling**:
- SMOTE-NC (Synthetic Minority Over-sampling with categorical features)
- Applied after CV split to prevent data leakage

**Model Architecture**:
- XGBoost with randomized hyperparameter search (50 iterations)
- Optimized for average precision (AUPRC)
- Tuned parameters: max_depth, colsample_bytree, learning_rate, regularization (alpha/lambda)

**Model Interpretability**:
- **SHAP TreeExplainer** on test sets for each fold
- Aggregates mean absolute SHAP values across out-of-fold (OOF) predictions
- Identifies top contributing features and their directionality

**Output Structure** (per assay):
```
output/TOX21_[ASSAY_NAME]/
├── run_summary.json              # CV metrics (ROC-AUC, PR-AUC, F1, etc.)
├── fold_details.csv              # Per-fold performance + class distributions
├── oof_predictions.csv           # Out-of-fold predictions for threshold tuning
├── feature_importance_summary.csv # XGBoost feature importances
├── feature_ranking.csv           # Ranked features by selection frequency
└── shap/
    ├── oof_shap.feather         # Sample-level SHAP values (OOF)
    ├── shap_global_rank_oof.csv # Global feature ranking by mean |SHAP|
    └── per_fold/                # Individual fold SHAP results
```

**Technologies**: Python 3.10, XGBoost, Boruta, SHAP, imbalanced-learn, scikit-learn

---

### 6. Results Analysis & Visualization (`06_results_analyses.ipynb`)

**Purpose**: Generate publication-quality figures and performance tables.

**Analyses**:

1. **Table 1: Model Performance Summary**
   - Cross-validated metrics for 41 assays (ROC-AUC, AUPRC, sensitivity, specificity, MCC)
   - Mean ± SD across 3 folds
   - Exported to `output/table_01_performance_table_paper.csv`

2. **Figure 1: Feature Reduction & Dataset Characteristics**
   - Waterfall plot showing Boruta feature selection stability
   - Heatmap of dataset properties (sample size, class balance, missingness)
   - Color-coded by cross-validation stability

3. **Figure 2: Model Performance Analysis**
   - Panel A: Assay-level AUPRC with error bars
   - Panel B: Fold-level AUPRC vs. sensitivity scatter plot
   - Panel C: Class balance visualization
   - Threshold line (AUPRC > 0.75) identifies 13 predictable assays

4. **Figure 3: SHAP Feature Importance Violin Plots**
   - Generated for all predictable assays (AUPRC > 0.75)
   - Top 20 features by mean absolute SHAP value
   - Color-coded by feature value (blue=low, red=high)
   - Side panel shows mean |SHAP| importance scores

**Technologies**: Python, matplotlib, seaborn, pandas

---

## Key Results

### Model Performance
- **13 of 41 assays** achieved AUPRC > 0.75 (considered "predictable")
- **Best performing assays**:
  - PXR Agonist: AUPRC 0.899 ± 0.020
  - PR Antagonist: AUPRC 0.892 ± 0.014
  - DT40 assays: AUPRC 0.872–0.889
- **High ROC-AUC across board**: Mean ROC-AUC > 0.85 for 30/41 assays

### Feature Insights
- **Transcriptomic features** consistently dominate importance rankings
- **Pathway activities** (PROGENy, MSigDB) more predictive than individual genes
- **Dose-response aggregation** (max, AUC_neg) captures non-linear chemical effects
- **Chemical fingerprints** provide complementary structural information

---

## Requirements

### Python Environment
```bash
# Core dependencies
python>=3.10
numpy>=1.26.4
pandas>=2.2.2
scikit-learn>=1.5.2
xgboost>=2.1.1
shap>=0.42.2
imbalanced-learn>=0.13.0
boruta>=0.4.3
rdkit>=2023.9.1
sqlalchemy>=2.0.40
feather-format>=0.4.1
joblib>=1.3.2
```

### R Environment
```bash
# R 4.3.3 with Bioconductor packages
DESeq2>=1.42.1
GSVA>=1.50.5
decoupleR>=2.8.0
progeny>=1.24.0
dorothea>=1.14.1
msigdbr>=25.1.1
BiocParallel>=1.36.0
tidyverse>=2.0.0
arrow>=21.0.0.1
```

---

## Data Availability

### Input Data Sources
1. **InvitroDB v4.2**: Tox21 HTS assay data (requires credentials)
2. **HTTr Count Data**:
   - MCF7: [GSE272548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272548)
   - U2OS: [Clowder](https://clowder.edap-cluster.com/files/6595a852e4b063812d5c77ef?dataset=61147fefe4b0856fdc65639b&space=&folder=65959df8e4b063812d5c7706)
   - HepRG: [GSE284321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284321)

### Intermediate Files (Not Included)
Due to file size constraints, the following intermediate files are not stored in this repository:
- `aggregated_results_[CELL_LINE].feather` (gene-level fold changes, ~500 MB each)
- Individual fold SHAP files (can be regenerated by running pipeline)

---

## Usage

### Running the Full Pipeline
```bash
# 1. Extract HTS assay data (requires InvitroDB access)
python analyses/01_assay_database_extraction_and_selection.py

# 2. Run differential expression analyses (R)
Rscript -e "rmarkdown::render('analyses/02a_MCF7_DE_analysis.Rmd')"
Rscript -e "rmarkdown::render('analyses/02b_U2OS_DE_analysis.Rmd')"
Rscript -e "rmarkdown::render('analyses/02c_HepRG_DE_analysis.Rmd')"

# 3. Combine and map transcriptomic activities (R)
Rscript -e "rmarkdown::render('analyses/03_combining_and_activity_mapping.Rmd')"

# 4. Integrate chemical structures and aggregate features (Jupyter)
jupyter notebook analyses/04_chemical_intergration_and_data_aggregation.ipynb

# 5. Train models with nested CV (Python, parallelized across 40 cores)
python analyses/05_modelling_and_evaluation.py

# 6. Generate figures and tables (Jupyter)
jupyter notebook analyses/06_results_analyses.ipynb
```

### Analyzing Individual Assays
```python
import pandas as pd
import json

# Load model results for a specific assay
assay_name = "TOX21_PXR_LUC_Agonist"
summary = json.load(open(f"output/{assay_name}/run_summary.json"))
shap_data = pd.read_feather(f"output/{assay_name}/shap/oof_shap.feather")

print(f"CV AUPRC: {summary['cv_pr_auc_mean']:.3f} ± {summary['cv_pr_auc_std']:.3f}")
```

---

## Citation

If you use this code or approach in your research, please cite:

```bibtex
@article{yourpaper2025,
  title={Mechanistically Interpretable Toxicity Prediction Through Multimodal Integration of Structure and Transcriptomics},
  author={Your Name et al.},
  journal={Journal Name},
  year={2025},
  note={GitHub: https://github.com/gcattebeke/toxicogenomics_paper}
}
```

---

## License

This project is licensed under the MIT License - see LICENSE file for details.

---

## Contact

For questions or collaboration inquiries, please open an issue on GitHub or contact the corresponding author.

---

## Acknowledgments

- **EPA CompTox Dashboard** for InvitroDB access and chemical annotations
- **GEO Database** for public HTTr data access
- **MSigDB, DoRothEA, PROGENy** developers for pathway analysis tools
- Computational resources provided by [Your Institution]

---

**Note**: Not all generated intermediate files are stored in this repository due to size constraints. The pipeline can regenerate these files from source data.

