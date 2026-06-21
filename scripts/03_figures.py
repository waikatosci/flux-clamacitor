#!/usr/bin/env python3
"""
03_figures.py
=============
Publication-quality figures for the WTP arsenic analysis (Figure 2 + ES&T
Supporting Information) of Hartland et al. (in review). Outputs both 600-dpi
PNG and vector PDF for each figure into figures/.

Depends on the intermediate result files produced by 02_wtp_arsenic_analysis.py
(run that first). Reads from data/, writes to figures/.

Figures
  S1  Treatment performance over time: inflow/outflow As + MAV, removal %.
  S2  Spearman correlation heatmap (raw-water drivers x removal/As metrics).
  S3  PCA biplot (PC1 vs PC2) coloured by As removal %, vectors = loadings.
  S4  Driver attribution: standardised OLS coefficients + RF permutation imp.
  S5  Mechanistic scatter: removal % vs Ca and vs Fe, coloured by era.

Style targets ACS / ES&T: Arial, single-column 3.33in or double 7.0in width,
minimal chartjunk, colour-blind-safe palette.

Run from the repository root:
    python scripts/02_wtp_arsenic_analysis.py   # produces intermediates
    python scripts/03_figures.py                # produces figures
"""
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data"          # repo convention: data + intermediates live in data/
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# ---------------------------------------------------------------- ES&T style
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 120, "savefig.dpi": 600, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02, "pdf.fonttype": 42, "ps.fonttype": 42,
})
# colour-blind-safe (Okabe-Ito)
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "red": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9",
      "yellow": "#F0E442", "grey": "#999999"}

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)
    print(f"  saved figures/{name}.png + .pdf")

df = pd.read_csv(OUT / "hcc_arsenic_tidy.csv", parse_dates=["datetime"])
res = json.load(open(OUT / "wtp_analysis_results.json"))
ann = pd.read_csv(OUT / "wtp_annual_summary.csv")
paired = df["As_tot_in"].notna() & df["As_tot_out"].notna()
dp = df[paired].copy()
dp = dp[~dp.get("qa_exclude_multivar", pd.Series(False, index=dp.index)).fillna(False)].copy()
dp["era"] = np.where(dp["coagulant"].astype(str).str.contains("ACH", na=False),
                     "Alum+ACH", "Alum")

print("Generating figures...")

# ================================================================ FIG S1
fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.8), sharex=True,
                         gridspec_kw={"height_ratios": [1.4, 1]})
ax = axes[0]
ax.scatter(dp["datetime"], dp["As_tot_in"]*1000, s=10, c=OI["blue"],
           alpha=0.55, edgecolors="none", label="Raw intake (total As)")
out_cens = df["As_tot_out_cens"] == True if "As_tot_out_cens" in df else pd.Series(False, index=df.index)
oc = df[paired & out_cens]; on = df[paired & ~out_cens]
ax.scatter(on["datetime"], on["As_tot_out"]*1000, s=10, c=OI["red"],
           alpha=0.6, edgecolors="none", label="Treated outlet (total As)")
ax.scatter(oc["datetime"], oc["As_tot_out"]*1000, s=10, facecolors="none",
           edgecolors=OI["red"], linewidths=0.6, alpha=0.6, label="Outlet < LOD")
ax.axhline(10, color="k", ls="--", lw=0.8)
ax.text(df["datetime"].min(), 10.4, "NZ MAV 0.010 mg L$^{-1}$", fontsize=6.5, va="bottom")
ax.set_ylabel("Total arsenic (µg L$^{-1}$)")
ax.set_yscale("log")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False,
          ncol=3, handletextpad=0.3, columnspacing=1.4, borderaxespad=0)
ax.text(0.01, 0.97, "a", transform=ax.transAxes, fontweight="bold", va="top")

ax2 = axes[1]
ax2.scatter(dp["datetime"], dp["As_removal_pct"], s=9, c=OI["grey"],
            alpha=0.5, edgecolors="none")
ax2.plot(pd.to_datetime(ann["year"], format="%Y"), ann["removal_med"],
         "-o", color=OI["green"], ms=3.5, lw=1.2, label="Annual median")
# Theil-Sen trend
s = dp[["datetime", "As_removal_pct"]].dropna().sort_values("datetime")
t = (s["datetime"] - s["datetime"].min()).dt.days / 365.25
sen = stats.theilslopes(s["As_removal_pct"].values, t.values)
ax2.plot(s["datetime"], sen[1] + sen[0]*t.values, color=OI["orange"], lw=1.3,
         label=f"Theil–Sen ({sen[0]:+.2f}%/yr, p<0.001)")
ax2.set_ylabel("As removal (%)")
ax2.set_xlabel("Year")
ax2.set_ylim(60, 100)
ax2.legend(loc="lower left", frameon=False)
ax2.text(0.01, 0.97, "b", transform=ax2.transAxes, fontweight="bold", va="top")
fig.align_ylabels(axes)
fig.subplots_adjust(top=0.90)
save(fig, "FigS1_performance_timeseries")

# ================================================================ FIG S2
metrics = ["As_removal_pct", "As_tot_in", "As_tot_out"]
drivers = ["pH", "temperature", "alkalinity", "Ca_tot", "Mg_tot", "Fe_tot",
           "Mn_tot", "DOC", "TOC", "turbidity", "hardness", "cond", "Cl",
           "SO4", "F", "algae", "Al_tot"]
drivers = [d for d in drivers if d in dp.columns]
allv = metrics + drivers
C = np.full((len(metrics), len(drivers)), np.nan)
P = np.full_like(C, np.nan)
for i, m in enumerate(metrics):
    for j, d in enumerate(drivers):
        sub = dp[[m, d]].dropna()
        if len(sub) >= 10:
            C[i, j], P[i, j] = stats.spearmanr(sub[m], sub[d])
fig, ax = plt.subplots(figsize=(7.0, 2.6))
im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(drivers)))
ax.set_xticklabels(drivers, rotation=45, ha="right")
ax.set_yticks(range(len(metrics)))
ax.set_yticklabels(["As removal %", "Inflow As", "Outflow As"])
for i in range(len(metrics)):
    for j in range(len(drivers)):
        if not np.isnan(C[i, j]):
            star = "*" if P[i, j] < 0.05 else ""
            ax.text(j, i, f"{C[i,j]:.2f}{star}", ha="center", va="center",
                    fontsize=5.6, color="white" if abs(C[i, j]) > 0.55 else "black")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("Spearman ρ", fontsize=7)
cb.ax.tick_params(labelsize=6)
ax.set_title("Rank correlation: As metrics vs raw-water drivers (* p<0.05)", fontsize=8)
save(fig, "FigS2_correlation_heatmap")

# ================================================================ FIG S3  PCA biplot
scores = pd.read_csv(OUT / "wtp_pca_scores.csv")
load = pd.read_csv(OUT / "wtp_pca_loadings.csv", index_col=0)
ev = res["pca"]["explained_var"]
fig, ax = plt.subplots(figsize=(4.6, 4.2))
sc = ax.scatter(scores["PC1"], scores["PC2"], c=scores["As_removal_pct"],
                cmap="viridis", s=18, alpha=0.8, edgecolors="none")
# scale arrows to span of scores for readability
span = np.percentile(np.abs(scores[["PC1", "PC2"]].values), 95)
scale = span / np.abs(load[["PC1", "PC2"]].values).max() * 0.9
# simple label de-overlap: sort by angle, nudge radially
for var in load.index:
    x, y = load.loc[var, "PC1"]*scale, load.loc[var, "PC2"]*scale
    ax.arrow(0, 0, x, y, color=OI["red"], width=0.003,
             head_width=0.06*span, length_includes_head=True, alpha=0.8, lw=0.5)
    ax.annotate(var, (x, y), xytext=(x*1.18, y*1.18), color=OI["red"],
                fontsize=6, ha="center", va="center")
ax.axhline(0, color="grey", lw=0.5, ls=":"); ax.axvline(0, color="grey", lw=0.5, ls=":")
lim = span * 1.5
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel(f"PC1 ({100*ev[0]:.0f}%)")
ax.set_ylabel(f"PC2 ({100*ev[1]:.0f}%)")
cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
cb.set_label("As removal (%)", fontsize=7); cb.ax.tick_params(labelsize=6)
ax.set_title(f"PCA of raw-water co-variates (n={res['pca']['n_samples']})", fontsize=8.5)
save(fig, "FigS3_PCA_biplot")

# ================================================================ FIG S4  attribution
# pretty variable labels (remove underscores / format subscripts)
LBL = {"As_tot_in": "As (intake)", "Ca_tot": "Ca", "Fe_tot": "Fe",
       "Mn_tot": "Mn", "alkalinity": "Alkalinity", "temperature": "Temperature",
       "pH": "pH", "DOC": "DOC"}
def pretty(v):
    return LBL.get(v, v.replace("_", " "))

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
# OLS standardised coefs (excl const), with significance
coefs = {k: v for k, v in res["ols"]["coefs"].items() if k != "const"}
pvs = res["ols"]["pvalues"]
ks = list(coefs.keys()); vals = [coefs[k] for k in ks]
cols = [OI["blue"] if pvs[k] < 0.05 else OI["grey"] for k in ks]
order = np.argsort(vals)
axes[0].barh([pretty(ks[i]) for i in order], [vals[i] for i in order],
             color=[cols[i] for i in order], edgecolor="k", linewidth=0.4)
axes[0].axvline(0, color="k", lw=0.7)
axes[0].set_xlabel("Standardised OLS coefficient")
axes[0].set_title(f"     Multiple regression (R²={res['ols']['r2']:.2f}, n={res['ols']['n']})",
                  fontsize=7.5, loc="left")
axes[0].text(0.0, 1.02, "A", transform=axes[0].transAxes, fontsize=7.5,
             fontweight="bold", va="bottom", ha="left")
axes[0].text(0.98, 0.04, "filled = p<0.05", transform=axes[0].transAxes,
             fontsize=6, ha="right", color=OI["grey"])
# RF importance
rf = pd.DataFrame(res["rf"]["permutation_importance"])
rf["label"] = rf["var"].map(pretty)
axes[1].barh(rf["label"][::-1], rf["importance"][::-1],
             xerr=rf["std"][::-1], color=OI["green"], edgecolor="k",
             linewidth=0.4, error_kw={"lw": 0.6})
axes[1].set_xlabel("Permutation importance (ΔR²)")
axes[1].set_title(f"     Random forest (n={res['rf']['n']})", fontsize=7.5, loc="left")
axes[1].text(0.0, 1.02, "B", transform=axes[1].transAxes, fontsize=7.5,
             fontweight="bold", va="bottom", ha="left")
save(fig, "FigS4_driver_attribution")

# ================================================================ FIG S5 mechanistic
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2))
era_col = {"Alum": OI["blue"], "Alum+ACH": OI["orange"]}
for ax, xvar, xlab in [(axes[0], "Ca_tot", "Total calcium (mg L$^{-1}$)"),
                       (axes[1], "Fe_tot", "Total iron (mg L$^{-1}$)")]:
    for era, g in dp.groupby("era"):
        sub = g[[xvar, "As_removal_pct"]].dropna()
        ax.scatter(sub[xvar], sub["As_removal_pct"], s=14, alpha=0.6,
                   c=era_col[era], edgecolors="none", label=era)
    sub = dp[[xvar, "As_removal_pct"]].dropna()
    if len(sub) > 10:
        rho, p = stats.spearmanr(sub[xvar], sub["As_removal_pct"])
        b = stats.theilslopes(sub["As_removal_pct"].values, sub[xvar].values)
        # draw over central 95% of x to avoid leverage from extreme tail
        xlo, xhi = np.percentile(sub[xvar], [1, 95])
        xs = np.linspace(xlo, xhi, 50)
        ax.plot(xs, b[1] + b[0]*xs, color="k", lw=1.1)
        ax.text(0.04, 0.06, f"ρ={rho:+.2f}, p<0.001" if p < 0.001 else
                f"ρ={rho:+.2f}, p={p:.3f}", transform=ax.transAxes, fontsize=6.5)
    ax.set_xlabel(xlab); ax.set_ylabel("As removal (%)")
    ax.set_ylim(60, 100)
axes[0].legend(frameon=False, loc="lower right", title="Coagulant era",
               title_fontsize=6.5)
axes[0].text(0.01, 0.99, "a", transform=axes[0].transAxes, fontweight="bold", va="top")
axes[1].text(0.01, 0.99, "b", transform=axes[1].transAxes, fontweight="bold", va="top")
save(fig, "FigS5_Ca_Fe_mechanism")

print("All figures complete.")
