#!/home/gcattebeke/miniconda3/envs/ML_env/bin/python
# ------------------------------ imports -------------------------------------- #
import contextlib
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import SMOTEN
from imblearn.over_sampling import SMOTENC
from imblearn.pipeline import Pipeline as ImbPipeline
from sklearn.model_selection import StratifiedGroupKFold

# ------------------------------ config --------------------------------------- #
DATA_PATH = Path("../../data/chemical_httr_assay_aggregated.feather")
OUT_BASE_ORIGINAL = Path("../../output")
OUT_BASE_REVISION = Path(".")
OUT_BASE_REVISION.mkdir(parents=True, exist_ok=True)

N_SPLITS = 3
RANDOM_STATE = 21
MIN_CV_PR_AUC = 0.75

# This revision script is intentionally baseline-only.
ABLATION_MODE = "baseline"

# If True, existing revised outputs are replaced.
OVERWRITE = False

# Cap SHAP background rows to reduce runtime/memory.
SHAP_BACKGROUND_MAX_ROWS = 10000

# Parallel settings.
TOTAL_CORES_BUDGET = 60
ENABLE_PARALLEL_ASSAYS = True
MAX_PARALLEL_ASSAYS = 13


def _safe_run_tag(selected_assay: str) -> str:
    assay_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", selected_assay)
    return assay_safe


def _parse_json_dict(value) -> dict:
    if isinstance(value, dict):
        return value
    if isinstance(value, str) and value.strip():
        try:
            parsed = json.loads(value)
            return parsed if isinstance(parsed, dict) else {}
        except json.JSONDecodeError:
            return {}
    return {}


def _parse_json_list(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, str) and value.strip():
        try:
            parsed = json.loads(value)
            return parsed if isinstance(parsed, list) else []
        except json.JSONDecodeError:
            return []
    return []


def _build_sampler(feature_cols: list[str], maccs_cols: list[str]):
    all_categorical = len(feature_cols) > 0 and all(c.startswith("MACCS_") for c in feature_cols)
    categorical_idx = [feature_cols.index(c) for c in feature_cols if c in maccs_cols]

    if all_categorical:
        return SMOTEN(sampling_strategy="auto", random_state=RANDOM_STATE)
    if categorical_idx:
        return SMOTENC(
            categorical_features=categorical_idx,
            sampling_strategy="auto",
            random_state=RANDOM_STATE,
        )
    return SMOTE(sampling_strategy="auto", random_state=RANDOM_STATE)


def _sample_background(
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    max_rows: int,
    seed: int,
) -> pd.DataFrame:
    if len(x_train) <= max_rows:
        return x_train.copy()

    rs = np.random.RandomState(seed)
    idx_pos = np.where(y_train == 1)[0]
    idx_neg = np.where(y_train == 0)[0]

    if len(idx_pos) == 0 or len(idx_neg) == 0:
        chosen = rs.choice(np.arange(len(x_train)), size=max_rows, replace=False)
        return x_train.iloc[chosen].copy()

    frac_pos = len(idx_pos) / len(y_train)
    n_pos = int(round(max_rows * frac_pos))
    n_pos = max(1, min(n_pos, len(idx_pos) - 1 if len(idx_pos) > 1 else 1))
    n_neg = max_rows - n_pos
    n_neg = max(1, min(n_neg, len(idx_neg) - 1 if len(idx_neg) > 1 else 1))

    # Ensure exact row count when class sizes are small.
    while n_pos + n_neg < max_rows:
        if n_pos < len(idx_pos):
            n_pos += 1
        elif n_neg < len(idx_neg):
            n_neg += 1
        else:
            break

    sel_pos = rs.choice(idx_pos, size=n_pos, replace=False)
    sel_neg = rs.choice(idx_neg, size=n_neg, replace=False)
    chosen = np.concatenate([sel_pos, sel_neg])
    rs.shuffle(chosen)

    return x_train.iloc[chosen].copy()


def _process_assay(
    selected_assay: str,
    assay_idx: int,
    total_assays: int,
    cv_pr_auc_mean: float,
    xgb_n_jobs: int,
    df_all: pd.DataFrame,
    httr_cols_all: list[str],
    chemical_cols_all: list[str],
) -> str:
    run_tag = _safe_run_tag(selected_assay)
    src_dir = OUT_BASE_ORIGINAL / run_tag
    src_fold_details = src_dir / "fold_details.csv"

    dst_dir = OUT_BASE_REVISION / run_tag
    dst_shap_per_fold = dst_dir / "shap" / "per_fold"
    dst_shap_per_fold.mkdir(parents=True, exist_ok=True)

    revision_marker = dst_dir / "revision_summary.json"
    if revision_marker.exists() and not OVERWRITE:
        return f"[{assay_idx}/{total_assays}] Skip {selected_assay}: already revised"

    log_path = dst_dir / "revision_log.txt"
    with open(log_path, "w") as log:
        httr_cols = list(httr_cols_all)
        chemical_cols = list(chemical_cols_all)

        df_full = df_all.copy()
        if "cell_type" in df_full.columns:
            df_full = df_full.drop(columns=["cell_type"])

        df_filtered = df_full.dropna(subset=[selected_assay]).copy()

        cols_to_drop = [
            c
            for c in df_filtered.columns
            if c.startswith("TOX21_") and c != selected_assay
        ]
        if "chnm" in df_filtered.columns:
            cols_to_drop.append("chnm")
        if cols_to_drop:
            df_filtered.drop(columns=cols_to_drop, inplace=True)

        ids_to_keep = [c for c in ["outcome_id", "epa_sample_id"] if c in df_filtered.columns]

        y_raw = df_filtered[selected_assay].astype(int).values
        groups = df_filtered["outcome_id"].astype(str).values

        cv = StratifiedGroupKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=RANDOM_STATE,
        )

        src_fold_df = pd.read_csv(src_fold_details)
        if "fold" not in src_fold_df.columns:
            return f"[{assay_idx}/{total_assays}] Skip {selected_assay}: fold_details.csv has no fold column"

        fold_rows = {int(row["fold"]): row for _, row in src_fold_df.iterrows()}
        per_fold_outputs = []
        successful_folds = 0

        for fold, (train_idx, test_idx) in enumerate(cv.split(df_filtered, y_raw, groups), start=1):
            if fold not in fold_rows:
                print(f"Fold {fold}: no row in fold_details.csv, skipping.", file=log)
                continue

            row = fold_rows[fold]

            train = df_filtered.iloc[train_idx].reset_index(drop=True)
            test = df_filtered.iloc[test_idx].reset_index(drop=True)

            y_train = train[selected_assay].astype(int).values
            y_test = test[selected_assay].astype(int).values

            if len(np.unique(y_train)) < 2:
                print(f"Fold {fold}: train is single class, skipping.", file=log)
                continue

            meta_cols = [c for c in ["outcome_id", "epa_sample_id"] if c in train.columns]
            x_train_df = train.drop(columns=[selected_assay] + meta_cols).copy()
            x_test_df = test.drop(columns=[selected_assay] + meta_cols).copy()

            maccs_cols = [c for c in x_train_df.columns if c.startswith("MACCS_")]
            if maccs_cols:
                x_train_df[maccs_cols] = x_train_df[maccs_cols].astype("int8")
                x_test_df[maccs_cols] = x_test_df[maccs_cols].astype("int8")

            selected_httr_features = _parse_json_list(row.get("selected_httr_features", "[]"))
            selected_httr_features = [
                c for c in selected_httr_features if c in x_train_df.columns and c in httr_cols
            ]

            feature_cols = chemical_cols + selected_httr_features
            feature_cols = [c for c in feature_cols if c in x_train_df.columns]

            if not feature_cols:
                print(f"Fold {fold}: no usable features, skipping.", file=log)
                continue

            x_train_selected = x_train_df[feature_cols]
            x_test_selected = x_test_df[feature_cols]

            best_params = _parse_json_dict(row.get("best_params", "{}"))

            model = xgb.XGBClassifier(
                objective="binary:logistic",
                eval_metric="aucpr",
                random_state=RANDOM_STATE,
                n_jobs=xgb_n_jobs,
                tree_method="hist",
                **best_params,
            )

            sampler = _build_sampler(feature_cols=feature_cols, maccs_cols=maccs_cols)
            pipeline = ImbPipeline([
                ("smote", sampler),
                ("classifier", model),
            ])

            pipeline.fit(x_train_selected, y_train)
            y_proba = pipeline.predict_proba(x_test_selected)[:, 1]

            xgb_model = pipeline.named_steps["classifier"]

            x_background = _sample_background(
                x_train=x_train_selected,
                y_train=y_train,
                max_rows=SHAP_BACKGROUND_MAX_ROWS,
                seed=RANDOM_STATE + fold,
            )

            try:
                explainer = shap.Explainer(
                    xgb_model,
                    x_background,
                    algorithm="tree",
                    feature_names=feature_cols,
                    model_output="probability",
                )

                with open(os.devnull, "w") as null_out, contextlib.redirect_stdout(null_out), contextlib.redirect_stderr(null_out):
                    shap_explanation = explainer(x_test_selected)

                shap_values_test = shap_explanation.values
                base_values_test = shap_explanation.base_values

                shap_cols = [f"SHAP_{c}" for c in feature_cols]
                shap_df_fold = pd.DataFrame(shap_values_test, columns=shap_cols)

                meta_df = test[ids_to_keep].copy() if ids_to_keep else pd.DataFrame(index=np.arange(len(test)))
                meta_df["fold"] = fold
                meta_df["y_true"] = y_test
                meta_df["y_proba"] = y_proba
                meta_df["base_value"] = base_values_test

                feat_vals_df = x_test_selected.reset_index(drop=True).copy()
                feat_vals_df.columns = feature_cols

                shap_out_fold = pd.concat(
                    [
                        meta_df.reset_index(drop=True),
                        feat_vals_df,
                        shap_df_fold,
                    ],
                    axis=1,
                )

                fold_path = dst_shap_per_fold / f"shap_test_fold_{fold}_unaug_bg.feather"
                shap_out_fold.to_feather(fold_path)

                mean_abs_shap = np.abs(shap_values_test).mean(axis=0)
                rank_fold = pd.DataFrame(
                    {
                        "fold": fold,
                        "feature": feature_cols,
                        "mean_abs_shap": mean_abs_shap,
                    }
                ).sort_values("mean_abs_shap", ascending=False)
                rank_path = dst_shap_per_fold / f"shap_rank_fold_{fold}_unaug_bg.csv"
                rank_fold.to_csv(rank_path, index=False)

                per_fold_outputs.append(shap_out_fold)
                successful_folds += 1

                print(
                    f"Fold {fold}: ok | n_features={len(feature_cols)} | bg_rows={len(x_background)}",
                    file=log,
                )

            except Exception as exc:
                print(f"Fold {fold}: SHAP failed: {exc}", file=log)
                continue

        if per_fold_outputs:
            shap_dir = dst_dir / "shap"
            shap_dir.mkdir(parents=True, exist_ok=True)

            oof_shap_df = pd.concat(per_fold_outputs, axis=0, ignore_index=True)
            oof_shap_path = shap_dir / "oof_shap_unaug_bg.feather"
            oof_shap_df.to_feather(oof_shap_path)

            shap_cols_all = [c for c in oof_shap_df.columns if c.startswith("SHAP_")]
            feature_names = [c.replace("SHAP_", "") for c in shap_cols_all]
            mean_abs = oof_shap_df[shap_cols_all].abs().mean(axis=0).values

            global_rank = pd.DataFrame(
                {
                    "feature": feature_names,
                    "mean_abs_shap_oof": mean_abs,
                }
            ).sort_values("mean_abs_shap_oof", ascending=False)

            global_rank["feature_type"] = global_rank["feature"].apply(
                lambda x: "httr"
                if any(x.endswith(suffix) for suffix in ["_max", "_dose_at_max", "_AUC_neg"])
                else "chemical"
                if x.startswith("MACCS_")
                else "other"
            )

            global_rank_path = shap_dir / "shap_global_rank_oof_unaug_bg.csv"
            global_rank.to_csv(global_rank_path, index=False)

        revision_summary = {
            "target_assay": selected_assay,
            "ablation_mode": ABLATION_MODE,
            "cell_type_filter": None,
            "selection_filter": f"cv_pr_auc_mean > {MIN_CV_PR_AUC}",
            "source_cv_pr_auc_mean": float(cv_pr_auc_mean),
            "source_run_dir": str(src_dir),
            "revision_run_dir": str(dst_dir),
            "n_splits": N_SPLITS,
            "random_state": RANDOM_STATE,
            "shap_background": "unaugmented_train_fold",
            "shap_background_max_rows": SHAP_BACKGROUND_MAX_ROWS,
            "xgb_n_jobs": xgb_n_jobs,
            "successful_folds": successful_folds,
            "overwrite": OVERWRITE,
        }
        with open(revision_marker, "w") as f:
            json.dump(revision_summary, f, indent=2)

    return f"[{assay_idx}/{total_assays}] Completed revision for {selected_assay}"


def main() -> None:
    df_all = pd.read_feather(DATA_PATH)
    assay_cols = [col for col in df_all.columns if col.startswith("TOX21_")]

    print(f"Found {len(assay_cols)} assays")
    print(f"Ablation mode: {ABLATION_MODE}")
    print(f"Filter: baseline cv_pr_auc_mean > {MIN_CV_PR_AUC}")
    print(f"Reading source outputs from: {OUT_BASE_ORIGINAL.resolve()}")
    print(f"Writing revised SHAP outputs to: {OUT_BASE_REVISION.resolve()}")

    httr_cols_all = [
        col
        for col in df_all.columns
        if col.endswith(("_max", "_dose_at_max", "_AUC_neg"))
    ]
    chemical_cols_all = [col for col in df_all.columns if col.startswith("MACCS_")]

    eligible = []
    for selected_assay in assay_cols:
        run_tag = _safe_run_tag(selected_assay)
        src_dir = OUT_BASE_ORIGINAL / run_tag
        src_run_summary = src_dir / "run_summary.json"
        src_fold_details = src_dir / "fold_details.csv"
        dst_dir = OUT_BASE_REVISION / run_tag
        revision_marker = dst_dir / "revision_summary.json"

        if revision_marker.exists() and not OVERWRITE:
            continue
        if not src_run_summary.exists() or not src_fold_details.exists():
            continue

        try:
            with open(src_run_summary, "r") as f:
                summary = json.load(f)
            cv_pr_auc_mean = summary.get("cv_pr_auc_mean", None)
            cv_pr_auc_mean = float(cv_pr_auc_mean) if cv_pr_auc_mean is not None else np.nan
        except Exception:
            cv_pr_auc_mean = np.nan

        if np.isfinite(cv_pr_auc_mean) and cv_pr_auc_mean > MIN_CV_PR_AUC:
            eligible.append((selected_assay, float(cv_pr_auc_mean)))

    print(f"Eligible baseline assays: {len(eligible)}")
    if not eligible:
        print("No eligible assays found. Done.")
        return

    if ENABLE_PARALLEL_ASSAYS and len(eligible) > 1:
        n_workers = min(MAX_PARALLEL_ASSAYS, len(eligible), TOTAL_CORES_BUDGET)
    else:
        n_workers = 1
    xgb_n_jobs = max(1, TOTAL_CORES_BUDGET // n_workers)

    print(
        f"Execution plan: workers={n_workers}, xgb_n_jobs_per_worker={xgb_n_jobs}, "
        f"total_core_budget={TOTAL_CORES_BUDGET}"
    )

    if n_workers == 1:
        for idx, (selected_assay, cv_pr_auc_mean) in enumerate(eligible, start=1):
            msg = _process_assay(
                selected_assay=selected_assay,
                assay_idx=idx,
                total_assays=len(eligible),
                cv_pr_auc_mean=cv_pr_auc_mean,
                xgb_n_jobs=xgb_n_jobs,
                df_all=df_all,
                httr_cols_all=httr_cols_all,
                chemical_cols_all=chemical_cols_all,
            )
            print(msg)
    else:
        futures = {}
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            for idx, (selected_assay, cv_pr_auc_mean) in enumerate(eligible, start=1):
                fut = executor.submit(
                    _process_assay,
                    selected_assay,
                    idx,
                    len(eligible),
                    cv_pr_auc_mean,
                    xgb_n_jobs,
                    df_all,
                    httr_cols_all,
                    chemical_cols_all,
                )
                futures[fut] = selected_assay

            for fut in as_completed(futures):
                selected_assay = futures[fut]
                try:
                    print(fut.result())
                except Exception as exc:
                    print(f"Failed {selected_assay}: {exc}")

    print("Done. Revised SHAP generation finished.")


if __name__ == "__main__":
    main()
