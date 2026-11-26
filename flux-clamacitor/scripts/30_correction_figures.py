"""
30_correction_figures.py
========================
Figures 3 and 4 and Table 1 of Hartland et al. (2026, ES&T) redrawn with the
values in the Addition and Correction:

  1. Ca at each site and date = median of the filtered and unfiltered ICP-MS
     replicates (four measurements), replacing the mean of the filtered pair.
  2. Shell bulk density = 0.632 g/cm3 (total dry mass / total ellipsoid volume,
     29 specimens, weighings corrected from the lab notes), replacing 1.15.

Everything else is as published: the published flux equations and flows,
positive-date averaging, alkalinity as measured, dv/dt = 0.0055 cm3/d and a
lakebed area of 7.41 km2.

The published figure scripts (04, 05) are run unchanged apart from their input
and output folders, so layout and styling match the paper.

Run from the repository root:
    python scripts/30_correction_figures.py
Writes figures/correction/Fig_3_mass_balance_corrected.{pdf,png},
figures/correction/Fig_4_density_corrected.{pdf,png} and
figures/correction/Table_1_corrected.csv
"""
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from clamacitor import io, fluxes, scaling                      # noqa: E402
from clamacitor.constants import LEGACY, Lake                  # noqa: E402
# (equivalent to the CORRECTION preset: LEGACY with median4 Ca and the corrected bulk density)

DATA = ROOT / "data"
WORK = ROOT / "output" / "correction_inputs"
FIG = ROOT / "figures" / "correction"
WORK.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SPEC = DATA / "correction" / "shell_mass_corrected_2026-11-26.xlsx"
ICPMS = DATA / "correction" / "icpms_replicates_2024-10_2025-07.csv"
PUBLISHED = DATA / "Full_data_for_Corbicula_Calcs_Landscape_Updated_16_1_26.xlsx"

STATION = {("Aniwaniwa Reserve", "Top"): "up", ("Golf Course", "Top"): "down",
           ("Lake Karapiro", "Top"): "top", ("Lake Karapiro", "Middle"): "mid",
           ("Lake Karapiro", "Bottom"): "bot"}
WIDE_SFX = {"up": "", "down": ".1", "bot": ".2", "mid": ".3", "top": ".4"}

# ------------------------------------------------------------------ inputs
rho = round(scaling.bulk_density_from_specimens(io.read_specimens(SPEC)), 3)
method = LEGACY.with_(bulk_density_t_m3=rho, ca_estimator="median4")
lake = Lake()
robust = io.ca_by_estimator(io.read_icpms_long(ICPMS), "median4")
robust_ug = robust.set_index(["date", "site"]).Ca_mgL * 1e3
print(f"bulk density {rho} g/cm3")

# long concentration sheet used by Figure 3, panels A and B
conc = pd.read_excel(DATA / "Full_data_for_Corbicula_Calcs.xlsx", sheet_name="Concs and Flow")
key = list(zip(pd.to_datetime(conc.Date).dt.normalize(),
               [STATION[(s, d)] for s, d in zip(conc.Station, conc.Depth)]))
new = robust_ug.reindex(key).values
assert np.isfinite(new).all(), "robust Ca missing for a published sample"
conc["Ca (ug/L)"] = new
with pd.ExcelWriter(WORK / "Full_data_for_Corbicula_Calcs.xlsx") as xw:
    conc.to_excel(xw, sheet_name="Concs and Flow", index=False)
shutil.copy(DATA / "Monthly_Water_Quality_Narrows_Summary.csv", WORK)

# mass-balance results table used by Figure 3 (E, F) and Figure 4
res = pd.read_csv(DATA / "Revised_Corbicula_Mass_Balance.csv")
res["Date"] = pd.to_datetime(res["Date"])
for site, sfx in WIDE_SFX.items():
    res[f"Ca (ug/L){sfx}"] = robust_ug.reindex(list(zip(res.Date, [site] * len(res)))).values

camp, flows = io.read_published_wide(PUBLISHED)
camp = camp.drop(columns="Ca_mgL").merge(robust, on=["date", "site"], how="left")
w = fluxes.to_wide(camp, flows).reindex(res.Date)
cols = {"S": ("longitudinal", fluxes.longitudinal), "R": ("vertical", fluxes.vertical)}
for tag, (_, fn) in cols.items():
    for a in ("Ca", "Alk"):
        J = fn(w, method, a) if tag == "R" else fn(w, method, a, lake)
        t = fluxes.t_caco3(J).values
        res[f"J_rem_{a}_{tag}"] = J.values
        res[f"m_CaCO3_{a}_{tag}"] = t
        res[f"V_bio_{a}_{tag}"] = scaling.biovolume(t, method)
        res[f"Density_{a}_{tag}"] = scaling.density(scaling.population(res[f"V_bio_{a}_{tag}"], method), lake)
res.to_csv(WORK / "Revised_Corbicula_Mass_Balance.csv", index=False)
res.to_csv(FIG / "mass_balance_corrected.csv", index=False)

# ------------------------------------------------------------------ Table 1
rows = []
pooled_N, pooled_D = [], []
for a in ("Ca", "Alk"):
    v = res.loc[res[f"V_bio_{a}_R"] > 0, f"V_bio_{a}_R"].values
    N = scaling.population(v, method)
    D = N / 7.41e6                       # area as used in the published table
    pooled_N.append(N); pooled_D.append(D)
    rows.append((a, N, D, 1))
N = np.concatenate(pooled_N); D = np.concatenate(pooled_D)
rows.append(("Ca + Alk", N, D, 0))       # published pooled spread used ddof = 0
tab = pd.DataFrame([{
    "method": "vertical depletion", "basis": b,
    "mean N (billions)": f"{Nn.mean() / 1e9:.2f} ± {Nn.std(ddof=dd) / 1e9:.2f}",
    "mean D (individual m-2)": f"{Dd.mean():.0f} ± {Dd.std(ddof=dd):.0f}",
    "median N (billions)": f"{np.median(Nn) / 1e9:.2f}",
    "median D (individual m-2)": f"{np.median(Dd):.0f}",
    "n dates": len(Nn)} for b, Nn, Dd, dd in rows])
tab.to_csv(FIG / "Table_1_corrected.csv", index=False)
print(tab.to_string(index=False))


# ------------------------------------------------------------------ figures
OFFSCALE = '''
YTOP_A, YTOP_B = 60, 1500


def mark_offscale(ax, groups, positions, top, fmt):
    """Mark values above the axis top with a small arrow and their value."""
    for vals, x in zip(groups, positions):
        for v in sorted(v for v in vals if v > top):
            ax.annotate(fmt.format(v), xy=(x, top), xytext=(x + 0.12, top * 0.93),
                        fontsize=PT_ANNOT if "PT_ANNOT" in globals() else 7, color="gray",
                        ha="left", va="center",
                        arrowprops=dict(arrowstyle="-|>", color="gray", lw=0.6))
'''
def run_published(script, out_stem, extra=()):
    src = (ROOT / "scripts" / script).read_text()
    src = src.replace('DATA = ROOT / "data"', f'DATA = Path(r"{WORK}")')
    src = src.replace('FIG = ROOT / "figures"', f'FIG = Path(r"{FIG}")')
    for a, b in extra:
        assert a in src, (script, a)
        src = src.replace(a, b)
    exec(compile(src, script, "exec"), {"__file__": str(ROOT / "scripts" / script),
                                        "__name__": "__main__"})


run_published("04_figure3_mass_balance.py", "Fig_3",
              extra=[('"Fig_3_mass_balance.pdf"', '"Fig_3_mass_balance_corrected.pdf"'),
                     ('"Fig_3_mass_balance.png"', '"Fig_3_mass_balance_corrected.png"')])
run_published("05_figure4_density.py", "Fig_4",
              extra=[('"Fig_4_density.pdf"', '"Fig_4_density_corrected.pdf"'),
                     ('"Fig_4_density.png"', '"Fig_4_density_corrected.png"'),
                     # the density scaling uses the whole lakebed, not the area < 1 m deep
                     ('"Res Ca\\n(<1m)", "Res Alk\\n(<1m)", "Combined\\nModel"]',
                      '"Res Ca", "Res Alk", "Combined\\nModel"]'),
                     ('["Res Ca\\n(<1m)", "Res Alk\\n(<1m)", "Combined\\nModel", "Survey\\nAvg"]',
                      '["Res Ca", "Res Alk", "Combined\\nModel", "Survey\\nAvg"]'),
                     # fixed axis tops; values above the top are marked at the edge with their value
                     ('ax.set_ylim(bottom=0)\nax.text(-0.09, 1.02, "A"',
                      'ax.set_ylim(0, YTOP_A)\nmark_offscale(ax, [v_bio_ca_r, v_bio_alk_r], [1, 2], YTOP_A, "{:.0f}")\nax.text(-0.09, 1.02, "A"'),
                     ('ax.set_ylim(bottom=0)\nax.text(-0.09, 1.02, "B"',
                      'ax.set_ylim(0, YTOP_B)\nmark_offscale(ax, [dens_ca_r, dens_alk_r, combined_densities], [1, 2, 3], YTOP_B, "{:,.0f}")\nax.text(-0.09, 1.02, "B"'),
                     ('# ------------------------------------------------------------------ 2. plotting',
                      OFFSCALE + '\n# ------------------------------------------------------------------ 2. plotting')])
