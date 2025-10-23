#!/usr/bin/env python3
"""
This script extracts assay data from InvitroDB v4.2 and applies a series of
filtering steps to prepare high-quality datasets for machine learning analysis.

Pipeline Steps:
    1. Database Extraction: Query InvitroDB for assay endpoints
    2. Pattern Filtering: Remove non-analytical file types
    3. Viability Filtering: Exclude cytotoxicity-driven hits
    4. Class Balance Filtering: Ensure sufficient samples per class
"""

import os                                                           # Python v3.10.16
import re                                                           # Python v3.10.16
import shutil                                                       # Python v3.10.16
from pathlib import Path                                            # Python v3.10.16
from typing import Iterable                                         # Python v3.10.16
import pandas as pd                                                 # v2.2.3                           
from sqlalchemy import create_engine, text, bindparam, Integer      # v2.0.40

# Configuration
BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = BASE_DIR / "output"
AEIDS_CSV = BASE_DIR / "input/TOX21_aeids_all.csv"

TEMP_DIR = BASE_DIR / ".temp_processing"
TEMP_00_DIR = TEMP_DIR / "00_raw_extracts"
TEMP_01_DIR = TEMP_DIR / "01_pattern_filtered"
TEMP_02_DIR = TEMP_DIR / "02_viability_filtered"
TEMP_03_DIR = TEMP_DIR / "03_class_balanced"

# Fit category codes (high-confidence labels per InvitroDB v4.2 documentation)
FITC_STRONG_ACTIVES: Iterable[int] = (41, 42, 37, 38)
FITC_STRONG_INACTIVES: Iterable[int] = (13, 15)

HITC_ACTIVE_THRESHOLD = 0.90
MIN_CLASS_COUNT = 200

EXCLUDE_PATTERNS = [
    r'ch[12]',      # Channel-specific data (only ratio files are used)
    r'Followup',    # Follow-up assays with sparse chemical coverage
    r'_RT_',        # Time-series data incompatible with endpoint analysis
    r'AutoFluor',   # Autofluorescence artifacts, not biological activity
]

DB_USER = os.getenv("DB_USER")
DB_PASS = os.getenv("DB_PASS")
DB_HOST = os.getenv("DB_HOST", "localhost")
DB_PORT = os.getenv("DB_PORT", "3306")
DB_NAME = os.getenv("DB_NAME", "invitrodb_v4_2")

engine = create_engine(
    f"mysql+mysqlconnector://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}",
    pool_pre_ping=True,
    future=True,
)

SQL_TEMPLATE = """
WITH filtered_mc5 AS (
    SELECT
        m5id,
        m4id,
        hitc,
        fitc
    FROM mc5
    WHERE hitc >= 0
      AND (
            fitc IN :fitc_actives
         OR fitc IN :fitc_inactives
      )
      AND aeid = :aeid
),
mc4_joined AS (
    SELECT
        f.m5id,
        f.m4id,
        f.hitc,
        f.fitc,
        mc4.spid
    FROM filtered_mc5 AS f
    INNER JOIN mc4 ON f.m4id = mc4.m4id
),
mc6_joined_without_flags AS (
    SELECT
        cj.m4id,
        cj.m5id,
        cj.spid,
        cj.hitc,
        cj.fitc
    FROM mc4_joined AS cj
    LEFT JOIN mc6 ON cj.m5id = mc6.m5id AND cj.m4id = mc6.m4id
    WHERE mc6.mc6_mthd_id IS NULL OR mc6.mc6_mthd_id = 20
),
annotated_results AS (
    SELECT 
        mw.hitc,
        mw.fitc, 
        s.spid,
        s.chid, 
        c.casn, 
        c.chnm,
        c.dsstox_substance_id AS dtxsid
    FROM mc6_joined_without_flags AS mw
    INNER JOIN sample   AS s ON mw.spid = s.spid
    INNER JOIN chemical AS c ON s.chid = c.chid
),
consensus AS (
    SELECT
        dtxsid,
        SUM(CASE WHEN hitc >= :thr THEN 1 ELSE 0 END) AS active_count,
        SUM(CASE WHEN hitc <  :thr THEN 1 ELSE 0 END) AS inactive_count
    FROM annotated_results
    GROUP BY dtxsid
)
SELECT
    dtxsid,
    CASE WHEN active_count >= inactive_count THEN 1 ELSE 0 END AS hitc
FROM consensus;
"""

def extract_database():
    """    
    Queries the database for each assay endpoint (aeid) and applies
    consensus hit-calling logic to resolve replicate measurements.
    """
    TEMP_00_DIR.mkdir(parents=True, exist_ok=True)

    aeids = pd.read_csv(AEIDS_CSV, header=0)

    stmt = text(SQL_TEMPLATE).bindparams(
        bindparam("fitc_actives", expanding=True, type_=Integer),
        bindparam("fitc_inactives", expanding=True, type_=Integer),
    )

    extracted = 0
    skipped = 0

    with engine.begin() as conn:
        for _, row in aeids.iterrows():
            aeid = int(row["aeid"])
            name = str(row["assay_component_endpoint_name"])
            safe_name = re.sub(r"[^A-Za-z0-9_\-]", "_", name)
            out_path = TEMP_00_DIR / f"{safe_name}.csv"

            if out_path.exists():
                skipped += 1
                continue

            params = {
                "aeid": aeid,
                "thr": float(HITC_ACTIVE_THRESHOLD),
                "fitc_actives": list(FITC_STRONG_ACTIVES),
                "fitc_inactives": list(FITC_STRONG_INACTIVES),
            }

            df = pd.read_sql(stmt, conn, params=params)
            df.to_csv(out_path, index=False)
            extracted += 1

    engine.dispose()
    print(f"Database extraction: {extracted} extracted, {skipped} skipped")


def filter_non_analytical_files():
    """    
    Removes files that are not suitable for endpoint analysis:
    """
    all_files = [f for f in os.listdir(TEMP_00_DIR) if f.endswith('.csv')]
    
    excluded = [
        f for f in all_files 
        if any(re.search(pattern, f, re.IGNORECASE) for pattern in EXCLUDE_PATTERNS)
    ]
    included = [f for f in all_files if f not in excluded]
    
    if TEMP_01_DIR.exists():
        shutil.rmtree(TEMP_01_DIR)
    TEMP_01_DIR.mkdir(parents=True, exist_ok=True)
    
    for f in included:
        shutil.copy(TEMP_00_DIR / f, TEMP_01_DIR / f)
    
    print(f"Pattern filtering: {len(excluded)}/{len(all_files)} files excluded")


def apply_viability_filter():
    """
    Filter out compounds flagged as toxic in viability assays.
    """
    viability_files = [f for f in os.listdir(TEMP_01_DIR) if 'viability' in f.lower()]
    all_files = [f for f in os.listdir(TEMP_01_DIR) if f.endswith('.csv')]
    
    TEMP_02_DIR.mkdir(parents=True, exist_ok=True)
    
    filtered_count = 0
    
    for viability_file in viability_files:
        try:
            viability_path = TEMP_01_DIR / viability_file
            viability_data = pd.read_csv(viability_path)
            
            if 'hitc' not in viability_data.columns or 'dtxsid' not in viability_data.columns:
                continue
            
            toxic_dtxsids = set(viability_data.loc[viability_data['hitc'] == 1, 'dtxsid'].tolist())
            
            lookup = viability_file.replace('_viability', '').replace('.csv', '')
            matching_files = [
                f for f in all_files 
                if lookup.lower() in f.lower() 
                and 'viability' not in f.lower()
                and f.lower().endswith('.csv')
            ]
            
            if not matching_files:
                continue

            for match_file in matching_files:
                match_path_01 = TEMP_01_DIR / match_file
                match_path_02 = TEMP_02_DIR / match_file
                match_data = pd.read_csv(match_path_01)
                
                if 'dtxsid' not in match_data.columns:
                    shutil.copy(match_path_01, match_path_02)
                    continue
                
                if toxic_dtxsids:
                    match_data_filtered = match_data[~match_data['dtxsid'].isin(toxic_dtxsids)]
                    match_data_filtered.to_csv(match_path_02, index=False)
                    filtered_count += 1
                else:
                    shutil.copy(match_path_01, match_path_02)
                    
        except Exception as e:
            print(f"Error processing viability file: {str(e)}")
    
    remaining_files = [f for f in all_files if 'viability' not in f.lower()]
    processed_files = set(os.listdir(TEMP_02_DIR))
    
    for f in remaining_files:
        if f not in processed_files:
            shutil.copy(TEMP_01_DIR / f, TEMP_02_DIR / f)
    
    print(f"Viability filtering: {filtered_count} assays processed")


def filter_low_count_files():
    """
    Filter out assays with insufficient class balance:
    """
    all_files = [f for f in os.listdir(TEMP_02_DIR) if f.endswith('.csv')]
    
    TEMP_03_DIR.mkdir(parents=True, exist_ok=True)
    
    excluded = []
    
    for file in all_files:
        df = pd.read_csv(TEMP_02_DIR / file)
        counts = df['hitc'].value_counts().to_dict()
        pos = counts.get(1, 0)
        neg = counts.get(0, 0)
        
        if pos + neg == 0 or pos < MIN_CLASS_COUNT or neg < MIN_CLASS_COUNT:
            excluded.append(file)
        else:
            shutil.copy(TEMP_02_DIR / file, TEMP_03_DIR / file)
    
    included_count = len(all_files) - len(excluded)
    
    print(f"Class balance filtering: {len(excluded)} assays excluded")
    print(f"Final dataset: {included_count} assays retained")


def merge_all_assays():
    """
    Merge all filtered assay files into a single CSV with dtxsid as rows
    and assay names as columns, with hitc values as cell values.
    """
    all_files = [f for f in os.listdir(TEMP_03_DIR) if f.endswith('.csv')]
    
    if not all_files:
        print("No assays to merge!")
        return
    
    merged_dfs = []
    
    for file in all_files:
        df = pd.read_csv(TEMP_03_DIR / file)
        assay_name = file.replace('.csv', '')
        df['assay_name'] = assay_name
        df = df[['dtxsid', 'hitc', 'assay_name']]
        merged_dfs.append(df)
    
    combined_df = pd.concat(merged_dfs, ignore_index=True)
    
    pivoted_df = combined_df.pivot_table(
        index='dtxsid',
        columns='assay_name',
        values='hitc',
        aggfunc='first'
    )
    
    pivoted_df.reset_index(inplace=True)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "all_assays_merged.csv"
    pivoted_df.to_csv(output_path, index=False)


def main():
    extract_database()
    filter_non_analytical_files()
    apply_viability_filter()
    filter_low_count_files()
    merge_all_assays()
    
    # cleanup
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR)

if __name__ == "__main__":
    main()
