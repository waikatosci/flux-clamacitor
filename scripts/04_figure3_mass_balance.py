#!/usr/bin/env python3
"""
04_figure3_mass_balance.py
==========================
Main-text Figure 3 of Hartland et al. (ES&T, es-2026-03421x): dissolved Ca and
alkalinity dynamics at Lake Karapiro, with CaCO3 removal fluxes from the
longitudinal and vertical (residence-time) depletion mass balances.

Replaces the notebook cell that wrote Fig_4_ECP_Final_Smoothed_Narrows_Final.png.
Two things changed for the ES&T production requirements:

  1. TRUE VECTOR OUTPUT.  The PDF is written by matplotlib's own PDF backend, so
     every line, patch and glyph stays as vector geometry.  Critically,
     pdf.fonttype = 42 embeds fonts as TrueType rather than matplotlib's default
     Type 3, which several publishers reject and which cannot be re-flowed or
     text-searched.  Do NOT produce the PDF by exporting a PNG and wrapping it
     (that is what the previously submitted Fig_3_mass_bal.pdf was: a single
     600-dpi raster inside a PDF container, with no live text at all).

  2. CANVAS SIZED TO THE JOURNAL COLUMN.  The notebook drew at 10 x 12 in
     (25.4 x 30.5 cm).  ACS sets figures to one column (~3.3 in) or two
     columns (~7.0 in), so a 25 cm-wide figure is reduced to roughly 72% on the
     page and every label shrinks with it.  Drawing at the final width instead
     means the point sizes below are the point sizes that get printed.

Run from the repository root:
    python scripts/04_figure3_mass_balance.py
Writes figures/Fig_3_mass_balance.pdf (vector, for submission) and .png
(600 dpi, for quick viewing / pasting into the manuscript file).
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline

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


CONC_FILE = data_file("Full_data_for_Corbicula_Calcs.xlsx")
RESULTS_FILE = data_file("Revised_Corbicula_Mass_Balance.csv")
SUMMARY_FILE = data_file("Monthly_Water_Quality_Narrows_Summary.csv")

# ------------------------------------------------------------------ page setup
# ACS: 1 column ~3.3 in, 2 columns ~7.0 in; figure + caption must fit the 9 in
# page depth, so 8.0 in of artwork leaves ~1 in for this figure's long caption.
FIG_W, FIG_H = 7.0, 8.0

# Point sizes below are FINAL PRINTED sizes because the canvas is drawn at the
# published width. ACS treats ~4.5 pt as the floor; these sit well clear of it.
PT_LABEL, PT_TICK, PT_LEGEND, PT_PANEL, PT_ANNOT = 9, 8, 8, 11, 7

mpl.rcParams.update({
    # DejaVu Sans ships with matplotlib, so the figure renders identically on
    # any machine, and it carries the subscript-three and sigma glyphs used in
    # the axis labels (some Helvetica builds do not).  Swap in Arial by putting
    # it first in the list; keep DejaVu Sans as the fallback.
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "mathtext.fontset": "dejavusans",
    "font.size": PT_LABEL,
    "axes.labelsize": PT_LABEL,
    "xtick.labelsize": PT_TICK,
    "ytick.labelsize": PT_TICK,
    "legend.fontsize": PT_LEGEND,
    "legend.title_fontsize": PT_LEGEND,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.5,
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

VLINE_DATE = pd.to_datetime("2025-05-01")
V_LAKE = 85e6  # m3, full lake volume

# ------------------------------------------- 1. historical stats, spline-smooth
df_stats = pd.read_csv(SUMMARY_FILE).set_index("Month").reindex(range(1, 13))
# Interpolate the numeric columns only; pandas >=3 raises on string columns
# (the Narrows summary carries a Month_Name text column).
_num = df_stats.select_dtypes(include=[np.number]).columns
df_stats = df_stats.copy()
df_stats[_num] = (df_stats[_num]
                  .interpolate(method="linear", limit_direction="both")
                  .ffill().bfill())

x_months = np.array([0] + list(range(1, 13)) + [13])


def _wrap(col):
    """Pad a 12-month series at both ends so the spline is periodic."""
    v = df_stats[col].values
    return np.concatenate(([df_stats.loc[12, col]], v, [df_stats.loc[1, col]]))


spl_ca_m = make_interp_spline(x_months, _wrap("Merged Ca_Mean"), k=3)
spl_ca_s = make_interp_spline(x_months, _wrap("Merged Ca_SD"), k=3)
spl_alk_m = make_interp_spline(x_months, _wrap("Alkalinity (mg/L) as CaCO3_Mean"), k=3)
spl_alk_s = make_interp_spline(x_months, _wrap("Alkalinity (mg/L) as CaCO3_SD"), k=3)


def get_smooth_hist(dates):
    d = pd.Series(pd.to_datetime(dates))
    m_dec = d.dt.month + (d.dt.day - 1) / d.dt.days_in_month
    return spl_ca_m(m_dec), spl_ca_s(m_dec), spl_alk_m(m_dec), spl_alk_s(m_dec)


# --------------------------------------------------- 2. observed concentrations
df_conc = pd.read_excel(CONC_FILE, sheet_name="Concs and Flow")
df_conc["DateTime"] = pd.to_datetime(df_conc["Date"])
df_conc["daily_date"] = pd.to_datetime(df_conc["DateTime"].dt.date)
df_conc["Ca_mg_L"] = df_conc["Ca (ug/L)"] / 1000

df_avg = (df_conc
          .groupby(["daily_date", "Station", "Depth"])
          .agg({"Ca_mg_L": "mean",
                "Alkalinity (mg/L)": "mean",
                "Flow (cumecs, m3/s)": "mean"})
          .reset_index())

ca_m_v, ca_s_v, alk_m_v, alk_s_v = get_smooth_hist(df_avg["daily_date"])
df_avg["Ca_Hist_Mean"], df_avg["Ca_Hist_SD"] = ca_m_v, ca_s_v
df_avg["Alk_Hist_Mean"], df_avg["Alk_Hist_SD"] = alk_m_v, alk_s_v

cols = ["daily_date", "Ca_mg_L", "Alkalinity (mg/L)",
        "Ca_Hist_Mean", "Ca_Hist_SD", "Alk_Hist_Mean", "Alk_Hist_SD"]

df_ani = df_avg.query("Station == 'Aniwaniwa Reserve' & Depth == 'Top'")[cols].sort_values("daily_date")
df_golf = df_avg.query("Station == 'Golf Course' & Depth == 'Top'")[cols].sort_values("daily_date")
df_ag = pd.merge(df_ani, df_golf, on="daily_date", how="inner", suffixes=("_ani", "_golf"))

df_top = df_avg.query("Station == 'Lake Karapiro' & Depth == 'Top'")[cols].sort_values("daily_date")
df_mid = df_avg.query("Station == 'Lake Karapiro' & Depth == 'Middle'")[cols].sort_values("daily_date")
df_bot = df_avg.query("Station == 'Lake Karapiro' & Depth == 'Bottom'")[cols].sort_values("daily_date")

df_k = pd.merge(df_top, df_bot, on="daily_date", how="inner", suffixes=("_top", "_bot"))
df_k = pd.merge(df_k, df_mid, on="daily_date", how="inner")
df_k = df_k.rename(columns={"Ca_mg_L": "Ca_mg_L_mid",
                            "Alkalinity (mg/L)": "Alkalinity_Mean_mid"})

# ---------------------------------------------------------- 3. mass-balance out
df_res = pd.read_csv(RESULTS_FILE)
df_res["Date"] = pd.to_datetime(df_res["Date"])
df_res["tau_days"] = V_LAKE / (df_res["Downstream Flow Karapiro (cumecs, m3/s)"] * 86400)


def add_season_labels(ax):
    starts = [(pd.Timestamp("2023-12-01"), "Summer"), (pd.Timestamp("2024-03-01"), "Autumn"),
              (pd.Timestamp("2024-06-01"), "Winter"), (pd.Timestamp("2024-09-01"), "Spring"),
              (pd.Timestamp("2024-12-01"), "Summer"), (pd.Timestamp("2025-03-01"), "Autumn"),
              (pd.Timestamp("2025-06-01"), "Winter"), (pd.Timestamp("2025-09-01"), "Spring")]
    trans = ax.get_xaxis_transform()
    x_min, x_max = ax.get_xlim()
    for i in range(len(starts) - 1):
        s_date, label = starts[i]
        e_date = starts[i + 1][0]
        s_num, e_num = mdates.date2num(s_date), mdates.date2num(e_date)
        mid = mdates.num2date(s_num + (e_num - s_num) / 2)
        if e_num > x_min and s_num < x_max:
            ax.text(mid, 1.02, label, transform=trans, ha="center", va="bottom",
                    fontsize=PT_ANNOT, clip_on=False)
            if s_num >= x_min:
                ax.text(s_date, 1.0, "|", transform=trans, ha="center", va="center",
                        fontsize=PT_ANNOT, clip_on=False)


def panel_letter(ax, letter):
    ax.text(-0.17, 1.04, letter, transform=ax.transAxes,
            fontsize=PT_PANEL, fontweight="bold", va="bottom")


# ------------------------------------------------------------------ 4. plotting
fig, axs = plt.subplots(3, 2, figsize=(FIG_W, FIG_H), sharex=True,
                        constrained_layout=True)

full_dates = pd.date_range(df_avg["daily_date"].min(), df_avg["daily_date"].max(), freq="D")
s_ca_m, s_ca_s, s_alk_m, s_alk_s = get_smooth_hist(full_dates)

# --- A: longitudinal Ca -------------------------------------------------------
ax = axs[0, 0]
ax.fill_between(full_dates, s_ca_m - s_ca_s, s_ca_m + s_ca_s, color="gray",
                alpha=0.15, label="Hist. Monthly ±1σ")
ax.plot(full_dates, s_ca_m, color="black", ls="--", lw=0.8, label="Hist. Median")
ax.plot(df_ag["daily_date"], df_ag["Ca_mg_L_ani"], label="ECP-S1", color="blue", marker="o")
ax.plot(df_ag["daily_date"], df_ag["Ca_mg_L_golf"], label="ECP-S3", color="green", marker="s")
ax.fill_between(df_ag["daily_date"], df_ag["Ca_mg_L_ani"], df_ag["Ca_mg_L_golf"],
                where=(df_ag["Ca_mg_L_ani"] > df_ag["Ca_mg_L_golf"]),
                color="blue", alpha=0.15, interpolate=True)
ax.fill_between(df_ag["daily_date"], df_ag["Ca_mg_L_ani"], df_ag["Ca_mg_L_golf"],
                where=(df_ag["Ca_mg_L_ani"] < df_ag["Ca_mg_L_golf"]),
                color="green", alpha=0.15, interpolate=True)
ax.axvline(VLINE_DATE, color="gray", lw=1)
ax.set_ylabel("Ca Concentration (mg L$^{-1}$)")
ax.legend(loc="upper right", framealpha=1.0, borderpad=0.2).get_frame().set_facecolor("white")
ax.set_ylim(top=8)
panel_letter(ax, "A")
add_season_labels(ax)

# --- B: depth-resolved Ca -----------------------------------------------------
ax = axs[0, 1]
ax.fill_between(full_dates, s_ca_m - s_ca_s, s_ca_m + s_ca_s, color="gray", alpha=0.15)
ax.plot(full_dates, s_ca_m, color="black", ls="--", lw=0.8)
ax.plot(df_k["daily_date"], df_k["Ca_mg_L_top"], label="2 m", color="blue", marker="o")
ax.plot(df_k["daily_date"], df_k["Ca_mg_L_mid"], label="13 m", color="darkgreen", marker="^")
ax.plot(df_k["daily_date"], df_k["Ca_mg_L_bot"], label="25 m", color="red", marker="s")
ax.fill_between(df_k["daily_date"], df_k["Ca_mg_L_top"], df_k["Ca_mg_L_bot"],
                where=(df_k["Ca_mg_L_top"] > df_k["Ca_mg_L_bot"]),
                color="red", alpha=0.15, interpolate=True)
ax.fill_between(df_k["daily_date"], df_k["Ca_mg_L_top"], df_k["Ca_mg_L_bot"],
                where=(df_k["Ca_mg_L_top"] < df_k["Ca_mg_L_bot"]),
                color="blue", alpha=0.15, interpolate=True)
ax.axvline(VLINE_DATE, color="gray", lw=1)
ax.set_ylabel("Ca Concentration (mg L$^{-1}$)")
ax.legend(loc="lower right", title="ECP-S2", framealpha=1.0,
          borderpad=0.2).get_frame().set_facecolor("white")
panel_letter(ax, "B")
add_season_labels(ax)

# --- C: longitudinal alkalinity ----------------------------------------------
ax = axs[1, 0]
ax.fill_between(full_dates, s_alk_m - s_alk_s, s_alk_m + s_alk_s, color="gray",
                alpha=0.15, label="Hist. Monthly ±1σ")
ax.plot(full_dates, s_alk_m, color="black", ls="--", lw=0.8)
ax.plot(df_ag["daily_date"], df_ag["Alkalinity (mg/L)_ani"], label="ECP-S1", color="blue", marker="o")
ax.plot(df_ag["daily_date"], df_ag["Alkalinity (mg/L)_golf"], label="ECP-S3", color="green", marker="s")
ax.fill_between(df_ag["daily_date"], df_ag["Alkalinity (mg/L)_ani"], df_ag["Alkalinity (mg/L)_golf"],
                where=(df_ag["Alkalinity (mg/L)_ani"] > df_ag["Alkalinity (mg/L)_golf"]),
                color="blue", alpha=0.15, interpolate=True)
ax.fill_between(df_ag["daily_date"], df_ag["Alkalinity (mg/L)_ani"], df_ag["Alkalinity (mg/L)_golf"],
                where=(df_ag["Alkalinity (mg/L)_ani"] < df_ag["Alkalinity (mg/L)_golf"]),
                color="green", alpha=0.15, interpolate=True)
ax.axvline(VLINE_DATE, color="gray", lw=1)
ax.set_ylabel("Alkalinity (mg L$^{-1}$ as CaCO₃)")
ax.legend(loc="lower left", framealpha=1.0, borderpad=0.2).get_frame().set_facecolor("white")
ax.set_ylim(30, 47.5)
panel_letter(ax, "C")

# --- D: depth-resolved alkalinity --------------------------------------------
ax = axs[1, 1]
ax.fill_between(full_dates, s_alk_m - s_alk_s, s_alk_m + s_alk_s, color="gray", alpha=0.15)
ax.plot(full_dates, s_alk_m, color="black", ls="--", lw=0.8)
ax.plot(df_k["daily_date"], df_k["Alkalinity (mg/L)_top"], label="2 m", color="blue", marker="o")
ax.plot(df_k["daily_date"], df_k["Alkalinity_Mean_mid"], label="13 m", color="darkgreen", marker="^")
ax.plot(df_k["daily_date"], df_k["Alkalinity (mg/L)_bot"], label="25 m", color="red", marker="s")
ax.fill_between(df_k["daily_date"], df_k["Alkalinity (mg/L)_top"], df_k["Alkalinity (mg/L)_bot"],
                where=(df_k["Alkalinity (mg/L)_top"] > df_k["Alkalinity (mg/L)_bot"]),
                color="red", alpha=0.15, interpolate=True)
ax.fill_between(df_k["daily_date"], df_k["Alkalinity (mg/L)_top"], df_k["Alkalinity (mg/L)_bot"],
                where=(df_k["Alkalinity (mg/L)_top"] < df_k["Alkalinity (mg/L)_bot"]),
                color="blue", alpha=0.15, interpolate=True)
ax.axvline(VLINE_DATE, color="gray", lw=1)
ax.text(VLINE_DATE - pd.Timedelta(days=10), 30.5, "Net CaCO₃ precipitation",
        rotation=90, va="bottom", ha="right", fontsize=PT_ANNOT,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1))
ax.text(VLINE_DATE + pd.Timedelta(days=10), 30.5, "Net CaCO₃ dissolution",
        rotation=90, va="bottom", ha="left", fontsize=PT_ANNOT,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1))
ax.set_ylabel("Alkalinity (mg L$^{-1}$ as CaCO₃)")
ax.set_ylim(30, 47.5)
panel_letter(ax, "D")

# --- E: longitudinal depletion fluxes ----------------------------------------
bar_w = pd.Timedelta("4D")
off_ca, off_alk = pd.Timedelta("-2D"), pd.Timedelta("2D")

ax = axs[2, 0]
ax.bar(df_res["Date"] + off_ca, df_res["m_CaCO3_Ca_S"], width=bar_w,
       label=r"Flux from $J_\mathrm{rem,\,Ca}$", color="tab:orange")
ax.bar(df_res["Date"] + off_alk, df_res["m_CaCO3_Alk_S"], width=bar_w, alpha=0.7,
       label=r"Flux from $J_\mathrm{rem,\,Alk}$", color="tab:purple")
ax.axvline(VLINE_DATE, color="gray", lw=1)
ax.set_ylabel(r"Removal Flux $J_{\mathrm{rem}}$, CaCO₃ (t d$^{-1}$)" + "\n" + "(Longitudinal Depletion)")
ax.legend(loc="upper right", framealpha=1.0, borderpad=0.2).get_frame().set_facecolor("white")
ax.set_ylim(bottom=0)
panel_letter(ax, "E")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

# --- F: vertical depletion fluxes + residence time ---------------------------
ax = axs[2, 1]
ax.bar(df_res["Date"] + off_ca, df_res["m_CaCO3_Ca_R"], width=bar_w,
       label=r"Flux from $J_\mathrm{rem,\,Ca}$", color="tab:orange")
ax.bar(df_res["Date"] + off_alk, df_res["m_CaCO3_Alk_R"], width=bar_w, alpha=0.7,
       label=r"Flux from $J_\mathrm{rem,\,Alk}$", color="tab:purple")
ax.axvline(VLINE_DATE, color="gray", lw=1)
ax.set_ylabel(r"Removal Flux $J_{\mathrm{rem}}$, CaCO₃ (t d$^{-1}$)" + "\n" + "(Vertical Depletion)")
ax2 = ax.twinx()
ax2.plot(df_res["Date"], df_res["tau_days"], color="black", ls="--", alpha=0.35,
         label="Residence Time")
ax2.set_ylabel(r"Residence Time $\tau$ (d)")
ax2.set_ylim(0, 10)
ax2.grid(False)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper right", framealpha=1.0,
          borderpad=0.2).get_frame().set_facecolor("white")
ax.set_ylim(bottom=0)
panel_letter(ax, "F")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

# ------------------------------------------------------------------- 5. output
pdf_path = FIG / "Fig_3_mass_balance.pdf"
fig.savefig(pdf_path, format="pdf",
            metadata={"Title": "Hartland et al. Figure 3 - Ca/alkalinity dynamics "
                               "and CaCO3 removal fluxes, Lake Karapiro",
                      "Creator": "flux-clamacitor 04_figure3_mass_balance.py"})
fig.savefig(FIG / "Fig_3_mass_balance.png", format="png", dpi=600)
plt.close(fig)

print(f"Wrote {pdf_path} ({FIG_W:.2f} x {FIG_H:.2f} in = "
      f"{FIG_W * 2.54:.1f} x {FIG_H * 2.54:.1f} cm), vector, TrueType-embedded")
print(f"Wrote {FIG / 'Fig_3_mass_balance.png'} (600 dpi companion)")
