#!/usr/bin/env python3
"""
05_figure4_density.py
=====================
Main-text Figure 4 of Hartland et al. (ES&T, es-2026-03421x): C. fluminea
biovolume formation rates from the residence-time mass balance (A) and the
population densities they imply, against the May 2025 benthic surveys (B).

Replaces the notebook cell that wrote Fig_4_density.png. Same two production
changes as 04_figure3_mass_balance.py:

  1. TRUE VECTOR OUTPUT via matplotlib's PDF backend, with pdf.fonttype = 42 so
     fonts embed as TrueType rather than matplotlib's default Type 3. The
     previously submitted Fig_4_density.pdf was a 600-dpi PNG wrapped in a PDF
     container - no live text, nothing production can rescale cleanly.

  2. CANVAS AT THE PUBLISHED WIDTH (two ACS columns, ~7.0 in) so the point sizes
     set below are the sizes that print. The notebook drew at 6 x 6 in and then
     cropped with bbox_inches='tight', which is why the supplied PDF came out at
     an arbitrary 15.53 cm rather than a column width. Space for the survey-site
     legend is reserved with subplots_adjust instead of tight-cropping, so the
     PDF page is exactly FIG_W x FIG_H.

Run from the repository root:
    python scripts/05_figure4_density.py
Writes figures/Fig_4_density.pdf (vector) and .png (600 dpi companion).
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def data_file(name):
    """Resolve a data file case-insensitively (repo filenames vary in case)."""
    exact = DATA / name
    if exact.exists():
        return exact
    for p in DATA.iterdir():
        if p.name.lower() == name.lower():
            return p
    raise FileNotFoundError(f"{name} not found in {DATA}")


RESULTS_FILE = data_file("Revised_Corbicula_Mass_Balance.csv")

# ------------------------------------------------------------------ page setup
FIG_W, FIG_H = 7.0, 6.2          # two ACS columns wide
LEGEND_RIGHT = 0.74              # fraction of width left for the axes
PT_LABEL, PT_TICK, PT_LEGEND, PT_PANEL = 9, 8, 8, 11

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "mathtext.fontset": "dejavusans",
    "font.size": PT_LABEL,
    "axes.labelsize": PT_LABEL,
    "axes.titlesize": PT_PANEL,
    "xtick.labelsize": PT_TICK,
    "ytick.labelsize": PT_TICK,
    "legend.fontsize": PT_LEGEND,
    "legend.title_fontsize": PT_LEGEND,
    "axes.linewidth": 0.8,
    "lines.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    # --- the settings that actually make the PDF acceptable to production ---
    "pdf.fonttype": 42,      # TrueType, not Type 3
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "pdf.compression": 6,
    "figure.dpi": 120,
    "savefig.dpi": 600,      # only affects the PNG companion
    "savefig.bbox": None,    # keep the PDF page exactly FIG_W x FIG_H
    "savefig.transparent": False,
})

# ----------------------------------------------------------- 1. load and derive
df_results = pd.read_csv(RESULTS_FILE)

v_bio_ca_r = df_results.loc[df_results["V_bio_Ca_R"] > 0, "V_bio_Ca_R"].tolist()
v_bio_alk_r = df_results.loc[df_results["V_bio_Alk_R"] > 0, "V_bio_Alk_R"].tolist()

DV_DT_MEDIAN_M3 = 0.0055 * 1e-6   # median shell growth, m3 d-1 individual-1
AREA_COLONIZABLE = 7.41e6         # m2


def compute_daily_density(bio_m3, dv_dt, area):
    """Daily density estimates (ind m-2) from biovolume formation rates."""
    return [(v / dv_dt) / area for v in bio_m3]


dens_ca_r = compute_daily_density(v_bio_ca_r, DV_DT_MEDIAN_M3, AREA_COLONIZABLE)
dens_alk_r = compute_daily_density(v_bio_alk_r, DV_DT_MEDIAN_M3, AREA_COLONIZABLE)
combined_densities = dens_ca_r + dens_alk_r

model_mean = float(np.mean(combined_densities))
model_std = float(np.std(combined_densities))

survey_data = {"Waipuke": 1134.5, "Horahora": 830.3, "Moana roa": 418.4,
               "Bobs landing": 394.7, "Little Waipa": 7.1}
survey_vals = list(survey_data.values())
survey_mean = float(np.mean(survey_vals))
survey_std = float(np.std(survey_vals, ddof=1))

print(f"--- Results summary (colonizable area {AREA_COLONIZABLE / 1e6:.3f} km2) ---")
print(f"Weighted model density: {model_mean:.0f} ± {model_std:.0f} ind m-2")
print(f"Survey mean density:    {survey_mean:.0f} ± {survey_std:.0f} ind m-2")

# ------------------------------------------------------------------ 2. plotting
fig, axs = plt.subplots(2, 1, figsize=(FIG_W, FIG_H))
fig.subplots_adjust(left=0.11, right=LEGEND_RIGHT, top=0.95, bottom=0.09, hspace=0.32)

box_kw = dict(widths=0.5,
              flierprops={"marker": "o", "markersize": 3,
                          "markerfacecolor": "gray", "markeredgecolor": "none",
                          "alpha": 0.5},
              medianprops={"color": "black"},
              boxprops={"linewidth": 0.8})

# --- A: biovolume formation rates --------------------------------------------
ax = axs[0]
ax.boxplot([v_bio_ca_r, v_bio_alk_r],
           tick_labels=["Residence\nCa", "Residence\nAlk"], **box_kw)
ax.set_ylabel(r"Biovolume Formation Rate (m$^3$ day$^{-1}$)")
ax.grid(True, linestyle="--", alpha=0.3)
ax.set_ylim(bottom=0)
ax.text(-0.09, 1.02, "A", transform=ax.transAxes,
        fontsize=PT_PANEL, fontweight="bold", va="bottom")

# --- B: implied densities vs surveys -----------------------------------------
ax = axs[1]
ax.boxplot([dens_ca_r, dens_alk_r, combined_densities],
           positions=[1, 2, 3],
           tick_labels=["Res Ca\n(<1m)", "Res Alk\n(<1m)", "Combined\nModel"],
           **box_kw)

# weighted model mean ±1 SD
pos_model = 3
ax.add_patch(mpl.patches.Rectangle((pos_model - 0.25, model_mean - model_std),
                                   0.5, 2 * model_std,
                                   color="blue", alpha=0.2, linewidth=0))
ax.plot([pos_model - 0.25, pos_model + 0.25], [model_mean, model_mean],
        color="blue", lw=1.6)

# survey mean ±1 SD in its own slot
pos_survey = 4
ax.add_patch(mpl.patches.Rectangle((pos_survey - 0.25, survey_mean - survey_std),
                                   0.5, 2 * survey_std,
                                   color="green", alpha=0.2, linewidth=0))
ax.plot([pos_survey - 0.25, pos_survey + 0.25], [survey_mean, survey_mean],
        color="green", lw=1.6)

ax.set_xlim(0.4, 4.6)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(["Res Ca\n(<1m)", "Res Alk\n(<1m)", "Combined\nModel", "Survey\nAvg"],
                   ha="center")
ax.set_ylabel(r"Estimated Density (ind m$^{-2}$)")
ax.grid(True, linestyle="--", alpha=0.3)
ax.set_ylim(bottom=0)
ax.text(-0.09, 1.02, "B", transform=ax.transAxes,
        fontsize=PT_PANEL, fontweight="bold", va="bottom")

site_colors = ["red", "blue", "green", "orange", "purple"]
for i, (site, val) in enumerate(survey_data.items()):
    ax.axhline(y=val, color=site_colors[i], linestyle="--", linewidth=0.8,
               alpha=0.7, label=site)
ax.legend(title="Survey Sites", loc="upper left", bbox_to_anchor=(1.02, 1),
          borderaxespad=0, frameon=False)

# ------------------------------------------------------------------- 3. output
pdf_path = FIG / "Fig_4_density.pdf"
fig.savefig(pdf_path, format="pdf",
            metadata={"Title": "Hartland et al. Figure 4 - C. fluminea biovolume "
                               "formation rates and implied population densities",
                      "Creator": "flux-clamacitor 05_figure4_density.py"})
fig.savefig(FIG / "Fig_4_density.png", format="png", dpi=600)
plt.close(fig)

print(f"Wrote {pdf_path} ({FIG_W:.2f} x {FIG_H:.2f} in = "
      f"{FIG_W * 2.54:.1f} x {FIG_H * 2.54:.1f} cm), vector, TrueType-embedded")
print(f"Wrote {FIG / 'Fig_4_density.png'} (600 dpi companion)")
