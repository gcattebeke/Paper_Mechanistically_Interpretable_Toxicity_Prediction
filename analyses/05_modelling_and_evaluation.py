#!/home/gcattebeke/miniconda3/envs/ML_env/bin/python

# ------------------------------ imports -------------------------------------- #
import numpy as np                                          # v1.26.4
import pandas as pd                                         # v2.2.2
from datetime import datetime                               # Python v3.10.13
from collections import Counter                             # Python v3.10.13
from sklearn.model_selection import StratifiedGroupKFold    # v1.5.2
from sklearn.metrics import (                               # v1.5.2
    roc_auc_score, precision_recall_curve, auc,
    accuracy_score, f1_score, confusion_matrix,
    average_precision_score, brier_score_loss, 
    log_loss)
from sklearn.ensemble import RandomForestClassifier         # v1.5.2
from boruta import BorutaPy                                 # v0.4.3     
import xgboost as xgb                                       # v2.1.1
from pathlib import Path                                    # Python v3.10.13 
import json                                                 # Python v3.10.13                                                    
import re                                                   # Python v3.10.13
from imblearn.over_sampling import SMOTENC                  # 0.13.0
from imblearn.over_sampling import SMOTE                    # 0.13.0
from imblearn.over_sampling import SMOTEN                   # 0.13.0
from sklearn.model_selection import RandomizedSearchCV      # v1.5.2
import shap                                                 # v0.42.2         
from imblearn.pipeline import Pipeline as ImbPipeline       # 0.13.0
import contextlib                                           # Python v3.10.13
import os                                                   # Python v3.10.13 
import sys                                                  # Python v3.10.13

# ------------------------------ config --------------------------------------- #
# Configuration
N_SPLITS = 3
RANDOM_STATE = 21
N_JOBS = 40
STABILITY_RUNS = 5
STABILITY_MIN_FRAC = 0.6

# ------------------------------ ablation settings ---------------------------- #

# options: "baseline", "structure_only", "httr_only", "single_cell_httr_structure", "double_cell_httr_structure"
ABLATION_MODE = "baseline"

# options: "MCF7", "U2OS", "HepRG"
SINGLE_CELL_TYPE = "U2OS" 
DOUBLE_CELL_TYPE_1 = "U2OS"  # options: "MCF7", "U2OS", "HepRG"
DOUBLE_CELL_TYPE_2 = "HepRG"  # options: "MCF7", "U2OS", "HepRG"

VALID_ABLATION_MODES = {
    "baseline",
    "structure_only",
    "httr_only",
    "single_cell_httr_structure",
    "double_cell_httr_structure"
}
if ABLATION_MODE not in VALID_ABLATION_MODES:
    raise ValueError(f"Unsupported ABLATION_MODE '{ABLATION_MODE}'. Allowed: {sorted(VALID_ABLATION_MODES)}")

if ABLATION_MODE == "baseline":
    OUT_BASE = Path("../output")
else:
    OUT_BASE = Path("../extra/ablation")
OUT_BASE.mkdir(parents=True, exist_ok=True)

# ------------------------------ Get all assays ------------------------------- #
chemical_httr_assay_aggregated = pd.read_feather("../data/chemical_httr_assay_aggregated.feather")
assay_cols = [col for col in chemical_httr_assay_aggregated.columns if col.startswith('TOX21_')]

print(f"\nFound {len(assay_cols)} assays to process\n")
print(f"Ablation mode: {ABLATION_MODE}")
if ABLATION_MODE == "single_cell_httr_structure":
    print(f"Cell type filter: {SINGLE_CELL_TYPE}")
elif ABLATION_MODE == "double_cell_httr_structure":
    print(f"Cell type filters: {DOUBLE_CELL_TYPE_1}, {DOUBLE_CELL_TYPE_2}")

# ------------------------------ Loop over all assays ------------------------- #
for assay_idx, selected_assay in enumerate(assay_cols, 1):
    ASSAY_SAFE = re.sub(r"[^A-Za-z0-9_.-]+", "_", selected_assay)
    run_tag = ASSAY_SAFE
    if ABLATION_MODE != "baseline":
        if ABLATION_MODE == "single_cell_httr_structure":
            run_tag = f"{ASSAY_SAFE}__{ABLATION_MODE}_{SINGLE_CELL_TYPE}"
        elif ABLATION_MODE == "double_cell_httr_structure":
            run_tag = f"{ASSAY_SAFE}__{ABLATION_MODE}_{DOUBLE_CELL_TYPE_1}_{DOUBLE_CELL_TYPE_2}"
        else:
            run_tag = f"{ASSAY_SAFE}__{ABLATION_MODE}"
    RUN_DIR = OUT_BASE / run_tag
    
    # Check if assay has already been processed
    if (RUN_DIR / 'run_summary.json').exists():
        print(f"--> Assay {selected_assay} already processed. Skipping...\n")
        continue
    
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    
    # Setup logging to file
    log_file = RUN_DIR / 'run_log.txt'
    log = open(log_file, 'w')
    
    print(f"[{assay_idx}/{len(assay_cols)}] Processing {selected_assay}...", flush=True)

    print(f"\n{'='*80}", file=log)
    print(f"Processing assay {assay_idx}/{len(assay_cols)}: {selected_assay}", file=log)
    print(f"{'='*80}\n", file=log)
    
    # Initialize storage for this assay
    fold_importances = []
    fold_metrics_detailed = []
    fold_selected_features = []
    fold_class_distributions = []
    shap_data_per_fold = []
    
    # ------------------------------ data ----------------------------------------- #
    httr_cols_all     = [col for col in chemical_httr_assay_aggregated.columns if col.endswith(('_max', '_dose_at_max', '_AUC_neg'))]
    chemical_cols_all = [col for col in chemical_httr_assay_aggregated.columns if col.startswith('MACCS_')]
    httr_cols = httr_cols_all
    chemical_cols = chemical_cols_all

    if ABLATION_MODE == "structure_only":
        httr_cols = []
    elif ABLATION_MODE == "httr_only":
        chemical_cols = []
    
    df_full = chemical_httr_assay_aggregated.copy()
    if ABLATION_MODE == "single_cell_httr_structure":
        if 'cell_type' not in df_full.columns:
            print(f"[{selected_assay}] Missing cell_type column, skipping.", file=log)
            log.close()
            continue
        df_full = df_full[df_full['cell_type'] == SINGLE_CELL_TYPE].copy()
    elif ABLATION_MODE == "double_cell_httr_structure":
        if 'cell_type' not in df_full.columns:
            print(f"[{selected_assay}] Missing cell_type column, skipping.", file=log)
            log.close()
            continue
        df_full = df_full[df_full['cell_type'].isin([DOUBLE_CELL_TYPE_1, DOUBLE_CELL_TYPE_2])].copy()

    if 'cell_type' in df_full.columns:
        df_full = df_full.drop(columns=['cell_type'])
    
    # Filter to samples with target assay data
    df_filtered = df_full.dropna(subset=[selected_assay]).copy()

    if df_filtered.empty:
        print(f"[{selected_assay}] No samples after ablation filtering, skipping.", file=log)
        log.close()
        continue
    
    # Keep only the selected assay and drop other assays
    cols_to_drop = [c for c in df_filtered.columns if c.startswith("TOX21_") and c != selected_assay]
    if 'chnm' in df_filtered.columns:
        cols_to_drop.append('chnm')
    df_filtered.drop(columns=cols_to_drop, inplace=True)
    
    IDs_TO_KEEP = [c for c in ['outcome_id', 'epa_sample_id'] if c in df_filtered.columns]
    
    y_raw = df_filtered[selected_assay].astype(int).values
    groups = df_filtered['outcome_id'].astype(str).values

    if len(np.unique(y_raw)) < 2:
        print(f"[{selected_assay}] Only one target class after filtering, skipping.", file=log)
        log.close()
        continue
    
    pos = int((y_raw == 1).sum())
    neg = int((y_raw == 0).sum())
    
    print(f"Target: {selected_assay}", file=log)
    print(f"Samples: {len(df_filtered):,} | Positives: {pos:,} | Negatives: {neg:,}\n", file=log)
    
    initial_stats = {
        'samples_with_target':  len(df_filtered),
        'n_httr_features':      len(httr_cols),
        'n_chemical_features':  len(chemical_cols),
        'target_positives':     pos,
        'target_negatives':     neg,
        'target_balance_ratio': min(pos, neg) / max(pos, neg) if max(pos, neg) > 0 else None,
        'unique_chemicals':     df_filtered['outcome_id'].nunique() if 'outcome_id' in df_filtered.columns else None,
        'ablation_mode':        ABLATION_MODE,
        'cell_type_filter':     SINGLE_CELL_TYPE if ABLATION_MODE == "single_cell_httr_structure" else None,
        'cell_type_filter_1':   DOUBLE_CELL_TYPE_1 if ABLATION_MODE == "double_cell_httr_structure" else None,
        'cell_type_filter_2':   DOUBLE_CELL_TYPE_2 if ABLATION_MODE == "double_cell_httr_structure" else None
    }
    
    # ------------------------------ CV setup ------------------------------------- #
    cv = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    
    feature_hits = Counter()
    shap_values_all_folds = []
    
    # OOF predictions
    oof_proba = np.full(len(df_filtered), np.nan)
    oof_true  = np.zeros(len(df_filtered), dtype=int)
    
    # ------------------------------ CV loop -------------------------------------- #
    for fold, (train_idx, test_idx) in enumerate(cv.split(df_filtered, y_raw, groups), start=1):
        train = df_filtered.iloc[train_idx].reset_index(drop=True)
        test  = df_filtered.iloc[test_idx].reset_index(drop=True)
        
        y_train = train[selected_assay].astype(int).values
        y_test  = test[selected_assay].astype(int).values
        
        # Safety check
        if len(np.unique(y_train)) < 2:
            print(f"[Fold {fold}] skipped (train has a single class).", file=log)
            continue
        
        # ----------------------- Prepare data ------------------------------------ #
        meta_cols = [c for c in ['outcome_id', 'epa_sample_id'] if c in train.columns]
        
        X_train_df = train.drop(columns=[selected_assay] + meta_cols).copy()
        X_test_df  = test.drop(columns=[selected_assay] + meta_cols).copy()
        
        maccs_cols = [c for c in X_train_df.columns if c.startswith("MACCS_")]
        
        if maccs_cols:
            X_train_df[maccs_cols] = X_train_df[maccs_cols].astype("int8")
            X_test_df[maccs_cols]  = X_test_df[maccs_cols].astype("int8")
        
        # ----------------------- Boruta feature selection (HTTr only) ------------ #
        httr_cols_train = [col for col in X_train_df.columns if col.endswith(('_max', '_dose_at_max', '_AUC_neg'))]
        httr_cols_train = [c for c in httr_cols_train if c in httr_cols]
        
        if httr_cols_train:
            rf_for_boruta = RandomForestClassifier(
                n_jobs=N_JOBS,
                n_estimators=1000,
                class_weight="balanced",
                random_state=RANDOM_STATE
            )
        
        if httr_cols_train:
            # stability-based feature selection over multiple random undersamples
            hits_counter = Counter()
            valid_runs   = 0
            pos_idx_all  = np.where(y_train == 1)[0]
            neg_idx_all  = np.where(y_train == 0)[0]

            for r in range(STABILITY_RUNS):
                if len(neg_idx_all) > len(pos_idx_all) and len(pos_idx_all) > 0:
                    rs = np.random.RandomState(RANDOM_STATE + r)
                    neg_sel = rs.choice(neg_idx_all, size=len(pos_idx_all), replace=False)
                    sel_idx = np.concatenate([pos_idx_all, neg_sel])
                else:
                    sel_idx = np.arange(len(y_train))
                
                X_train_httr_boruta = X_train_df.iloc[sel_idx][httr_cols_train].values
                y_train_boruta      = y_train[sel_idx]
                
                boruta = BorutaPy(
                    estimator=rf_for_boruta,
                    n_estimators='auto',
                    perc=100,
                    verbose=0,
                    random_state=RANDOM_STATE + r
                )
                
                try:
                    boruta.fit(X_train_httr_boruta, y_train_boruta)
                    for f, keep in zip(httr_cols_train, boruta.support_):
                        if keep:
                            hits_counter[f] += 1
                    valid_runs += 1
                except Exception:
                    continue
            
            # Keep features selected in >= STABILITY_MIN_FRAC of valid runs
            threshold_hits = int(np.ceil(STABILITY_MIN_FRAC * max(valid_runs, 1)))
            selected_httr_features = [f for f, h in hits_counter.items() if h >= threshold_hits]
            
            # Fallback if stability keeps nothing
            if len(selected_httr_features) == 0:
                rf_for_boruta.fit(X_train_df[httr_cols_train].values, y_train)
                importances = rf_for_boruta.feature_importances_
                order = np.argsort(importances)[::-1]
                top_k = max(25, min(200, int(np.sqrt(len(httr_cols_train)))))
                selected_httr_features = [httr_cols_train[i] for i in order[:top_k]]
        else:
            selected_httr_features = []
        
        feature_hits.update(selected_httr_features)
        fold_selected_features.append(selected_httr_features)
        
        # ----------------------- Prepare pipeline with SMOTE + XGB -------------- #
        feature_cols = chemical_cols + selected_httr_features
        if not feature_cols:
            print(f"[Fold {fold}] no features available after ablation filter.", file=log)
            continue

        # If features are all MACCS (categorical), use SMOTEN instead of SMOTE-NC.
        all_categorical = len(feature_cols) > 0 and all(c.startswith("MACCS_") for c in feature_cols)
        
        X_train_selected = X_train_df[feature_cols]
        X_test_selected  = X_test_df[feature_cols]
        
        categorical_idx = [feature_cols.index(c) for c in feature_cols if c in maccs_cols]
        
        # Store class distribution info before SMOTE
        pre_smote_pos = int((y_train == 1).sum())
        pre_smote_neg = int((y_train == 0).sum())
        
        base_model = xgb.XGBClassifier(
            objective='binary:logistic',
            eval_metric='aucpr',
            random_state=RANDOM_STATE,
            n_jobs=1,
            tree_method='hist'
        )
        
        # Pipeline with SMOTE + XGB
        if all_categorical:
            sampler = SMOTEN(
                sampling_strategy='auto',
                random_state=RANDOM_STATE
            )
        elif categorical_idx:
            sampler = SMOTENC(
                categorical_features=categorical_idx,
                sampling_strategy='auto',
                random_state=RANDOM_STATE
            )
        else:
            sampler = SMOTE(
                sampling_strategy='auto',
                random_state=RANDOM_STATE
            )

        pipeline = ImbPipeline([
            ('smote', sampler),
            ('classifier', base_model)
        ])
        
        param_distributions = {
            'classifier__max_depth': [3, 4, 5, 6, 7, 8],
            'classifier__colsample_bytree': [0.4, 0.6, 0.8, 1.0],
            'classifier__reg_alpha': [0, 1e-3, 1e-2, 1e-1, 1, 5],
            'classifier__reg_lambda': [0.5, 1, 2, 5, 10],
            'classifier__learning_rate': [0.02, 0.05, 0.1],
            'classifier__n_estimators': [400, 600, 800, 1000, 1200]
        }
        
        # Use StratifiedGroupKFold for inner CV
        train_groups = train['outcome_id'].astype(str).values
        inner_cv = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)
        
        tuner = RandomizedSearchCV(
            estimator=pipeline,
            param_distributions=param_distributions,
            n_iter=50,
            scoring='average_precision',
            cv=inner_cv,
            verbose=0,
            n_jobs=N_JOBS,
            refit=True,
            random_state=RANDOM_STATE
        )
        
        tuner.fit(X_train_selected, y_train, groups=train_groups)
        model = tuner.best_estimator_
        
        # Get SMOTE-resampled data info for class distributions
        X_resampled_check = model.named_steps['smote'].fit_resample(X_train_selected, y_train)
        post_smote_pos = int((X_resampled_check[1] == 1).sum())
        post_smote_neg = int((X_resampled_check[1] == 0).sum())
        
        # Store class distribution info 
        fold_class_distributions.append({
            'fold': fold,
            'train_samples': len(y_train),
            'test_samples': len(y_test),
            'train_pos_pre_smote': pre_smote_pos,
            'train_neg_pre_smote': pre_smote_neg,
            'train_pos_post_smote': post_smote_pos,
            'train_neg_post_smote': post_smote_neg,
            'test_pos': int((y_test == 1).sum()),
            'test_neg': int((y_test == 0).sum()),
            'smote_applied': True
        })
        
        # ----------------------- Evaluate ---------------------------------------- #
        y_proba   = model.predict_proba(X_test_selected)[:, 1]
        xgb_model = model.named_steps['classifier']
        
        # Store OOF predictions for global threshold selection
        oof_proba[test_idx] = y_proba
        oof_true[test_idx]  = y_test
        
        # ----------------------- SHAP values (save per-fold + aggregate) --------- #
        shap_dir = RUN_DIR / 'shap'
        shap_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            X_resampled    = model.named_steps['smote'].fit_resample(X_train_selected, y_train)[0]
            X_resampled_df = pd.DataFrame(X_resampled, columns=feature_cols)
            
            explainer = shap.Explainer(
                xgb_model, 
                X_resampled_df, 
                algorithm="tree", 
                feature_names=feature_cols,
                model_output="probability"
            )

            with open(os.devnull, "w") as f, contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                shap_explanation = explainer(X_test_selected)
            shap_values_test = shap_explanation.values
            base_values_test = shap_explanation.base_values
            
            # Store mean absolute SHAP values for aggregation
            mean_abs_shap = np.abs(shap_values_test).mean(axis=0)
            shap_values_all_folds.append(pd.DataFrame({
                'fold': fold,
                'feature': feature_cols,
                'mean_abs_shap': mean_abs_shap
            }))
            
            # Build per-fold SHAP dataframe with metadata, features, and SHAP values
            shap_cols = [f"SHAP_{c}" for c in feature_cols]
            shap_df_fold = pd.DataFrame(shap_values_test, columns=shap_cols)
            meta_df = test[IDs_TO_KEEP].copy() if IDs_TO_KEEP else pd.DataFrame(index=np.arange(len(test)))
            meta_df['fold'] = fold
            meta_df['y_true'] = y_test
            meta_df['y_proba'] = y_proba
            meta_df['base_value'] = base_values_test
            feat_vals_df = X_test_selected.reset_index(drop=True).copy()
            feat_vals_df.columns = feature_cols
            
            shap_out_fold = pd.concat([meta_df.reset_index(drop=True), 
                                      feat_vals_df,
                                      shap_df_fold], axis=1)
            shap_data_per_fold.append(shap_out_fold)
            
            # Save per-fold SHAP file
            (shap_dir / 'per_fold').mkdir(exist_ok=True, parents=True)
            shap_out_path = shap_dir / 'per_fold' / f'shap_test_fold_{fold}.feather'
            shap_out_fold.to_feather(shap_out_path)
            
            # Save per-fold feature ranking
            shap_rank_fold = pd.DataFrame({
                'fold': fold,
                'feature': feature_cols,
                'mean_abs_shap': mean_abs_shap
            }).sort_values('mean_abs_shap', ascending=False)
            shap_rank_path = shap_dir / 'per_fold' / f'shap_rank_fold_{fold}.csv'
            shap_rank_fold.to_csv(shap_rank_path, index=False)
            
        except Exception as e:
            print(f"[Fold {fold}] SHAP calculation failed: {e}", file=log)
        
        # ----------------------- Metrics ----------------------------------------- #
        try:
            roc = roc_auc_score(y_test, y_proba)
        except ValueError:
            roc = np.nan
        
        precision, recall, thr = precision_recall_curve(y_test, y_proba)
        pr_auc = auc(recall, precision) if len(recall) > 1 else np.nan
        
        try:
            ap = average_precision_score(y_test, y_proba)
        except ValueError:
            ap = np.nan
        
        brier = brier_score_loss(y_test, y_proba)
        
        try:
            ll = log_loss(y_test, y_proba)
        except ValueError:
            ll = np.nan
        
        # Per-fold optimal F1 threshold
        if thr.size > 0:
            f1s = (2 * precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-12)
            best_idx = int(np.nanargmax(f1s))
            fold_opt_threshold = float(thr[best_idx])
            
            y_pred_opt = (y_proba >= fold_opt_threshold).astype(int)
            acc_opt = accuracy_score(y_test, y_pred_opt)
            f1_opt = f1_score(y_test, y_pred_opt, zero_division=0)
            cm_opt = confusion_matrix(y_test, y_pred_opt, labels=[0, 1])
        else:
            fold_opt_threshold = 0.5
            y_pred_opt = (y_proba >= 0.5).astype(int)
            acc_opt = np.nan
            f1_opt = np.nan
            cm_opt = confusion_matrix(y_test, y_pred_opt, labels=[0, 1])
        
        # Store fold metrics
        fold_metrics_detailed.append({
            'fold': fold,
            'n_train': len(y_train),
            'n_test': len(y_test),
            'train_pos': int((y_train == 1).sum()),
            'train_neg': int((y_train == 0).sum()),
            'test_pos': int((y_test == 1).sum()),
            'test_neg': int((y_test == 0).sum()),
            'roc_auc': roc,
            'pr_auc': pr_auc,
            'average_precision': ap,
            'brier_score': brier,
            'log_loss': ll,
            'optimal_threshold': fold_opt_threshold,
            'accuracy_at_optimal': acc_opt,
            'f1_at_optimal': f1_opt,
            'tn': int(cm_opt[0, 0]),
            'fp': int(cm_opt[0, 1]),
            'fn': int(cm_opt[1, 0]),
            'tp': int(cm_opt[1, 1]),
            'precision_at_optimal': int(cm_opt[1, 1]) / (int(cm_opt[1, 1]) + int(cm_opt[0, 1])) if (int(cm_opt[1, 1]) + int(cm_opt[0, 1])) > 0 else 0,
            'recall_at_optimal': int(cm_opt[1, 1]) / (int(cm_opt[1, 1]) + int(cm_opt[1, 0])) if (int(cm_opt[1, 1]) + int(cm_opt[1, 0])) > 0 else 0,
            'specificity_at_optimal': int(cm_opt[0, 0]) / (int(cm_opt[0, 0]) + int(cm_opt[0, 1])) if (int(cm_opt[0, 0]) + int(cm_opt[0, 1])) > 0 else 0,
            'httr_features_selected': len(selected_httr_features),
            'best_params': json.dumps({k.replace('classifier__', ''): v for k, v in tuner.best_params_.items() if k.startswith('classifier__')})
        })
        
        # Store feature importances
        if hasattr(xgb_model, 'feature_importances_'):
            imp_df = pd.DataFrame({
                'fold': fold,
                'feature': feature_cols,
                'importance': xgb_model.feature_importances_
            })
            fold_importances.append(imp_df)
        
        # Print fold summary to log
        cmopt_str = f"[{cm_opt.ravel()[0]},{cm_opt.ravel()[1]};{cm_opt.ravel()[2]},{cm_opt.ravel()[3]}]"
        print(
            f"[Fold {fold}/{N_SPLITS}] "
            f"Boruta: {len(selected_httr_features):>3} | "
            f"ROC: {roc:.3f} | PR-AUC: {pr_auc:.3f} | AP: {ap:.3f} | "
            f"F1: {f1_opt:.3f} | CM: {cmopt_str}",
            file=log
        )

    # ------------------------------ Global threshold (OOF) ----------------------- #
    mask = ~np.isnan(oof_proba)
    y_true_oof = oof_true[mask]
    y_prob_oof = oof_proba[mask]

    if len(y_true_oof) == 0:
        global_threshold = 0.5
        print("\nNo valid OOF predictions found; using default threshold 0.5", file=log)
    else:
        prec_oof, rec_oof, thr_oof = precision_recall_curve(y_true_oof, y_prob_oof)
        if len(thr_oof) == 0:
            global_threshold = 0.5
            print("\nUnable to compute PR thresholds; using default threshold 0.5", file=log)
        else:
            f1s = (2 * prec_oof[:-1] * rec_oof[:-1]) / (prec_oof[:-1] + rec_oof[:-1] + 1e-12)
            best_idx = int(np.nanargmax(f1s))
            global_threshold = float(thr_oof[best_idx])
        
        y_pred_oof = (y_prob_oof >= global_threshold).astype(int)
        oof_roc = roc_auc_score(y_true_oof, y_prob_oof) if len(np.unique(y_true_oof)) > 1 else float('nan')
        oof_pr = auc(rec_oof, prec_oof) if len(rec_oof) > 1 else float('nan')
        oof_acc = accuracy_score(y_true_oof, y_pred_oof)
        oof_f1 = f1_score(y_true_oof, y_pred_oof, zero_division=0)
        oof_cm = confusion_matrix(y_true_oof, y_pred_oof, labels=[0, 1])
        
        print(f"\nOOF Results:", file=log)
        print(f"ROC-AUC: {oof_roc:.3f} | PR-AUC: {oof_pr:.3f} | "
              f"Accuracy: {oof_acc:.3f} | F1: {oof_f1:.3f}", file=log)
        print(f"Confusion Matrix: {oof_cm.tolist()}", file=log)
        print(f"Global Threshold: {global_threshold:.4f}\n", file=log)

    # ------------------------------ SHAP aggregation across folds ---------------- #
    try:
        per_fold_dir = RUN_DIR / 'shap' / 'per_fold'
        if per_fold_dir.exists():
            shap_files = sorted(per_fold_dir.glob('shap_test_fold_*.feather'))
            if shap_files:
                # Aggregate all per-fold SHAP data into OOF SHAP file
                oof_shap_df = pd.concat([pd.read_feather(p) for p in shap_files], axis=0, ignore_index=True)
                oof_shap_path = RUN_DIR / 'shap' / 'oof_shap.feather'
                oof_shap_df.to_feather(oof_shap_path)
                
                # Calculate global feature ranking from OOF SHAP
                shap_cols_all = [c for c in oof_shap_df.columns if c.startswith('SHAP_')]
                feature_names_from_shap = [c.replace('SHAP_', '') for c in shap_cols_all]
                mean_abs = oof_shap_df[shap_cols_all].abs().mean(axis=0).values
                global_rank_df = pd.DataFrame({
                    'feature': feature_names_from_shap,
                    'mean_abs_shap_oof': mean_abs
                }).sort_values('mean_abs_shap_oof', ascending=False)
                
                # Merge with Boruta hit counts
                try:
                    feat_hits_df = pd.DataFrame.from_dict(feature_hits, orient='index').reset_index()
                    feat_hits_df.columns = ['feature', 'hit_count']
                    global_rank_df = global_rank_df.merge(feat_hits_df, on='feature', how='left')
                except Exception:
                    pass
                
                # Add feature type classification
                global_rank_df['feature_type'] = global_rank_df['feature'].apply(
                    lambda x: 'httr' if any(x.endswith(suffix) for suffix in ['_max', '_dose_at_max', '_AUC_neg'])
                    else 'chemical' if x.startswith('MACCS_') else 'other'
                )
                
                # Save global ranking files
                (RUN_DIR / 'shap').mkdir(exist_ok=True, parents=True)
                global_rank_df.to_csv(RUN_DIR / 'shap' / 'shap_global_rank_oof.csv', index=False)
                
                # Aggregate all per-fold rankings
                rank_files = sorted(per_fold_dir.glob('shap_rank_fold_*.csv'))
                if rank_files:
                    all_ranks_df = pd.concat([pd.read_csv(p) for p in rank_files], axis=0, ignore_index=True)
                    all_ranks_df.to_csv(RUN_DIR / 'shap' / 'shap_rank_all_folds.csv', index=False)
    except Exception as _e:
        pass

    # ------------------------------ Save essential outputs ----------------------- #

    # 1. Run summary (comprehensive metadata and performance metrics)
    if fold_metrics_detailed:
        metrics_df = pd.DataFrame(fold_metrics_detailed)
        cv_summary = {
            'target_assay': selected_assay,
            'n_folds': N_SPLITS,
            'random_state': RANDOM_STATE,
            'global_threshold': global_threshold,
            'cv_roc_auc_mean': float(metrics_df['roc_auc'].mean()),
            'cv_roc_auc_std': float(metrics_df['roc_auc'].std()),
            'cv_pr_auc_mean': float(metrics_df['pr_auc'].mean()),
            'cv_pr_auc_std': float(metrics_df['pr_auc'].std()),
            'cv_average_precision_mean': float(metrics_df['average_precision'].mean()),
            'cv_average_precision_std': float(metrics_df['average_precision'].std()),
            'cv_f1_mean': float(metrics_df['f1_at_optimal'].mean()),
            'cv_f1_std': float(metrics_df['f1_at_optimal'].std()),
            'cv_accuracy_mean': float(metrics_df['accuracy_at_optimal'].mean()),
            'cv_accuracy_std': float(metrics_df['accuracy_at_optimal'].std()),
            'cv_precision_mean': float(metrics_df['precision_at_optimal'].mean()),
            'cv_precision_std': float(metrics_df['precision_at_optimal'].std()),
            'cv_recall_mean': float(metrics_df['recall_at_optimal'].mean()),
            'cv_recall_std': float(metrics_df['recall_at_optimal'].std()),
            'cv_specificity_mean': float(metrics_df['specificity_at_optimal'].mean()),
            'cv_specificity_std': float(metrics_df['specificity_at_optimal'].std()),
            'cv_brier_score_mean': float(metrics_df['brier_score'].mean()),
            'cv_brier_score_std': float(metrics_df['brier_score'].std()),
            'cv_log_loss_mean': float(metrics_df['log_loss'].mean()),
            'cv_log_loss_std': float(metrics_df['log_loss'].std()),
            'oof_roc_auc': float(oof_roc) if 'oof_roc' in locals() and not np.isnan(oof_roc) else None,
            'oof_pr_auc': float(oof_pr) if 'oof_pr' in locals() and not np.isnan(oof_pr) else None,
            'oof_accuracy': float(oof_acc) if 'oof_acc' in locals() else None,
            'oof_f1': float(oof_f1) if 'oof_f1' in locals() else None,
            **initial_stats,
        }
        
        with open(RUN_DIR / 'run_summary.json', 'w') as f:
            json.dump(cv_summary, f, indent=2)

    # 2. Consolidated fold details (metrics + class distributions + selected features)
    if fold_metrics_detailed and fold_class_distributions:
        fold_details = pd.DataFrame(fold_metrics_detailed)
        
        # Merge with class distributions
        class_dist_df = pd.DataFrame(fold_class_distributions)
        fold_details = fold_details.merge(
            class_dist_df[['fold', 'train_pos_pre_smote', 'train_neg_pre_smote', 
                           'train_pos_post_smote', 'train_neg_post_smote', 'smote_applied']],
            on='fold',
            how='left'
        )
        
        # Add selected features as a column
        fold_details['selected_httr_features'] = [
            json.dumps(feats) for feats in fold_selected_features
        ]
        
        fold_details.to_csv(RUN_DIR / 'fold_details.csv', index=False)

    # 3. Feature importance summary (aggregated across folds for main figures)
    if fold_importances:
        importances_df = pd.concat(fold_importances, axis=0, ignore_index=True)
        
        importance_stats = importances_df.groupby('feature')['importance'].agg([
            'count', 'mean', 'std', 'min', 'max'
        ]).reset_index()
        importance_stats.columns = ['feature', 'fold_count', 'importance_mean', 
                                    'importance_std', 'importance_min', 'importance_max']
        importance_stats = importance_stats.sort_values('importance_mean', ascending=False)
        
        # Add feature type classification
        importance_stats['feature_type'] = importance_stats['feature'].apply(
            lambda x: 'httr' if any(x.endswith(suffix) for suffix in ['_max', '_dose_at_max', '_AUC_neg'])
                     else 'chemical' if x.startswith('MACCS_')
                     else 'other'
        )
        
        # Add Boruta selection frequency
        feat_hits_df = pd.DataFrame.from_dict(feature_hits, orient='index', columns=['selection_count']).reset_index()
        feat_hits_df.columns = ['feature', 'selection_count']
        importance_stats = importance_stats.merge(feat_hits_df, on='feature', how='left')
        importance_stats['selection_count'] = importance_stats['selection_count'].fillna(0).astype(int)
        
        importance_stats.to_csv(RUN_DIR / 'feature_importance_summary.csv', index=False)

    # 4. SHAP analysis summary (fold-level statistics)
    if shap_values_all_folds:
        shap_df = pd.concat(shap_values_all_folds, axis=0, ignore_index=True)
        
        shap_summary = shap_df.groupby('feature')['mean_abs_shap'].agg([
            'count', 'mean', 'std'
        ]).reset_index()
        shap_summary.columns = ['feature', 'fold_count', 'mean_abs_shap', 'std_abs_shap']
        shap_summary = shap_summary.sort_values('mean_abs_shap', ascending=False)
        
        # Add feature type classification
        shap_summary['feature_type'] = shap_summary['feature'].apply(
            lambda x: 'httr' if any(x.endswith(suffix) for suffix in ['_max', '_dose_at_max', '_AUC_neg'])
                     else 'chemical' if x.startswith('MACCS_')
                     else 'other'
        )
        
        # Note: This is fold-level aggregation; sample-level OOF SHAP is in shap/oof_shap.feather
        shap_summary.to_csv(RUN_DIR / 'shap' / 'shap_summary_by_fold.csv', index=False)

    # 5. OOF predictions (basic predictions without SHAP)
    if len(y_true_oof) > 0:
        oof_df = df_filtered[IDs_TO_KEEP].copy() if IDs_TO_KEEP else pd.DataFrame(index=np.arange(len(df_filtered)))
        oof_df['y_true'] = oof_true
        oof_df['y_proba'] = oof_proba
        oof_df['y_pred'] = (oof_proba >= global_threshold).astype(int)
        
        # Add fold information
        fold_assignment = np.full(len(df_filtered), -1, dtype=int)
        for fold, (_, test_idx) in enumerate(cv.split(df_filtered, y_raw, groups), start=1):
            fold_assignment[test_idx] = fold
        oof_df['fold'] = fold_assignment
        
        oof_df = oof_df[oof_df['fold'] != -1]  # Remove any samples not in test sets
        oof_df.to_csv(RUN_DIR / 'oof_predictions.csv', index=False)

    # Close log file
    log.close()
    print(f"[{assay_idx}/{len(assay_cols)}] ✓ Completed {selected_assay}", flush=True)

print("\n" + "="*80)
print("All assays processed successfully!")
print("="*80)
