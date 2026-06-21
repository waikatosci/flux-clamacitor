#!/usr/bin/env python3
"""
02_wtp_arsenic_analysis.py
==========================
Driver-attribution analysis of the Hamilton (Waiora) WTP operational record,
supporting Figure 2 and the Supporting Information of:

  Hartland et al. (in review) "Carbonate-system perturbation by invasive
  bivalves threatens drinking-water security in calcium-limited catchments",
  Environmental Science & Technology.

Consumes the cleaned, tidy operational dataset (data/hcc_arsenic_tidy.csv;
see data/wtp_arsenic_qaqc_log.txt and data/wtp_arsenic_data_dictionary.csv)
and produces the intermediate result files consumed by 03_figures.py:

  data/wtp_annual_summary.csv     annual median As-removal % (+ counts)
  data/wtp_pca_scores.csv         PCA scores (PC1, PC2) + As_removal_pct
  data/wtp_pca_loadings.csv       PCA loadings per covariate
  data/wtp_analysis_results.json  PCA variance, OLS coefs/p/R2, RF importances

Methods (manuscript Section 2.7): Spearman rank correlations, PCA of the
standardised covariate matrix, OLS multiple regression (standardised
coefficients), and random-forest regression with permutation importance.

Run from the repository root:
    python scripts/02_wtp_arsenic_analysis.py
"""
import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
import statsmodels.api as sm

DATA = "data"
RANDOM_STATE = 42

# --- Multivariate sample-construction (recovered from the original pipeline) -
# The manuscript driver models (Fig. 2; R2 = 0.33, OLS reported n = 95) were
# built with a coverage-based covariate/row selection followed by median
# imputation and standardisation, NOT naive listwise deletion over all analytes:
#   1. Start from paired intake-outlet events with a measured target, minus the
#      qa_exclude_multivar flag.
#   2. Keep only covariates with >= COV_MIN_COVERAGE non-null coverage.
#   3. Keep only rows with >= ROW_MIN_PRESENT of those covariates present.
#   4. Median-impute remaining gaps, then z-standardise.
# This recovers the manuscript's covariate set (Ca, Fe, Mn, Al, TOC, algae, pH,
# temperature, alkalinity) and R2 ~ 0.33 with Ca the sole significant positive
# predictor and alkalinity/pH/temperature significant negatives.
COV_MIN_COVERAGE = 0.40   # covariate retained if >= this fraction non-null
ROW_MIN_PRESENT = 0.70    # row retained if >= this fraction of kept covs present

# ------------------------------------------------------------------ load
df = pd.read_csv(os.path.join(DATA, "hcc_arsenic_tidy.csv"),
                 parse_dates=["datetime"])

# Paired intake-outlet events only; honour the multivariate QA exclusion flag.
paired = df["As_tot_in"].notna() & df["As_tot_out"].notna()
dp = df[paired].copy()
if "qa_exclude_multivar" in dp.columns:
    dp = dp[~dp["qa_exclude_multivar"].fillna(False)].copy()
print(f"Paired events: {paired.sum()}  (after multivariate QA: {len(dp)})")

TARGET = "As_removal_pct"
DRIVERS = ["pH", "temperature", "alkalinity", "Ca_tot", "Mg_tot", "Fe_tot",
           "Mn_tot", "DOC", "TOC", "turbidity", "hardness", "cond", "Cl",
           "SO4", "F", "algae", "Al_tot"]
DRIVERS = [d for d in DRIVERS if d in dp.columns]

# ------------------------------------------------------------------ annual summary
ann = (dp.dropna(subset=[TARGET])
         .groupby("year")[TARGET]
         .agg(removal_med="median", removal_mean="mean", n="count")
         .reset_index())
ann.to_csv(os.path.join(DATA, "wtp_annual_summary.csv"), index=False)
print(f"  wrote data/wtp_annual_summary.csv ({len(ann)} years)")

# ------------------------------------------------------------------ build modelling matrix
# Coverage-based covariate/row selection (recovered original methodology; see
# header), then median imputation + standardisation.
base = dp.dropna(subset=[TARGET]).copy()
coverage = {c: base[c].notna().mean() for c in DRIVERS}
kept = [c for c in DRIVERS if coverage[c] >= COV_MIN_COVERAGE]
row_frac = base[kept].notna().mean(axis=1)
M = base[row_frac >= ROW_MIN_PRESENT].copy()

print(f"  covariates >= {COV_MIN_COVERAGE:.0%} coverage ({len(kept)}): {kept}")
print(f"  rows with >= {ROW_MIN_PRESENT:.0%} of those present: n = {len(M)}")

DRIVERS = kept   # downstream PCA / OLS / RF operate on the retained covariates
Xfilled = M[DRIVERS].fillna(M[DRIVERS].median(numeric_only=True))
X = Xfilled.values
y = M[TARGET].values
Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
n = len(M)

# ------------------------------------------------------------------ PCA
pca = PCA(random_state=RANDOM_STATE)
scores = pca.fit_transform(Xz)
ev = pca.explained_variance_ratio_

pc_scores = pd.DataFrame(scores[:, :2], columns=["PC1", "PC2"])
pc_scores[TARGET] = y
pc_scores.to_csv(os.path.join(DATA, "wtp_pca_scores.csv"), index=False)

loadings = pd.DataFrame(pca.components_[:2].T, index=DRIVERS,
                        columns=["PC1", "PC2"])
loadings.to_csv(os.path.join(DATA, "wtp_pca_loadings.csv"))
print("  wrote data/wtp_pca_scores.csv + data/wtp_pca_loadings.csv")

# ------------------------------------------------------------------ OLS (standardised)
yz = (y - y.mean()) / y.std(ddof=0)
Xc = sm.add_constant(pd.DataFrame(Xz, columns=DRIVERS))
ols = sm.OLS(yz, Xc).fit()
ols_res = {
    "coefs": ols.params.to_dict(),
    "pvalues": ols.pvalues.to_dict(),
    "r2": float(ols.rsquared),
    "n": int(n),
}

# ------------------------------------------------------------------ Random forest
rf = RandomForestRegressor(n_estimators=500, random_state=RANDOM_STATE,
                           n_jobs=-1)
rf.fit(Xz, y)
perm = permutation_importance(rf, Xz, y, n_repeats=50,
                              random_state=RANDOM_STATE, n_jobs=-1)
order = np.argsort(perm.importances_mean)[::-1]
rf_res = {
    "n": int(n),
    "oob_r2": None,
    "permutation_importance": [
        {"var": DRIVERS[i],
         "importance": float(perm.importances_mean[i]),
         "std": float(perm.importances_std[i])}
        for i in order
    ],
}

# ------------------------------------------------------------------ assemble JSON
results = {
    "n_paired_events": int(paired.sum()),
    "target": TARGET,
    "drivers": DRIVERS,
    "pca": {
        "n_samples": int(n),
        "explained_var": [float(v) for v in ev],
    },
    "ols": ols_res,
    "rf": rf_res,
}
with open(os.path.join(DATA, "wtp_analysis_results.json"), "w") as f:
    json.dump(results, f, indent=2)
print("  wrote data/wtp_analysis_results.json")
print(f"OLS R2 = {ols_res['r2']:.2f}; "
      f"Ca std-beta = {ols_res['coefs'].get('Ca_tot', float('nan')):+.2f}")
print("Analysis complete.")
