"""Input readers. Every reader returns the same long 'campaign' table:

    date, site, Ca_mgL, Alk_mgL

with site one of the codes in SITES. Flows are a separate daily table:

    date, Q_out_m3s, Q_in_m3s, level_m
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

# canonical site codes. 'up'/'down' are the longitudinal pair, 'top'/'bot'
# the vertical pair at the dam face.
SITES = {
    "up": "Aniwaniwa (S1)",
    "down": "Golf Course (S3)",
    "bot": "Karapiro dam, bottom",
    "mid": "Karapiro dam, middle",
    "top": "Karapiro dam, top",
    "narrows": "Narrows (Karapiro outflow)",
}
_ROBUST_SITE = {"Aniwaniwa": "up", "Golf_Course": "down", "Karapiro_bot": "bot",
                "Karapiro_mid": "mid", "Karapiro_top": "top", "Narrows": "narrows"}

# position of each station block in the published wide workbook
_WIDE_ORDER = ["up", "down", "bot", "mid", "top"]


# ---------------------------------------------------------------------------
# Published (legacy) workbook
# ---------------------------------------------------------------------------
def read_published_wide(path):
    """Read 'Full_data_for_Corbicula_Calcs_Landscape_Updated_16_1_26.xlsx'.

    Returns (campaign_long, flows) where flows holds the published Arapuni
    (Q_in) and Karapiro (Q_out) values for each sampling date."""
    w = pd.read_excel(path)
    w.columns = w.columns.str.strip()
    w["Date"] = pd.to_datetime(w["Date"])
    rows = []
    for i, site in enumerate(_WIDE_ORDER):
        sfx = "" if i == 0 else f".{i}"
        rows.append(pd.DataFrame({
            "date": w["Date"], "site": site,
            "Ca_mgL": w[f"Ca (ug/L){sfx}"] / 1e3,
            "Alk_mgL": w[f"Alkalinity (mg/L){sfx}"],
        }))
    camp = pd.concat(rows, ignore_index=True)
    flows = pd.DataFrame({
        "date": w["Date"],
        "Q_in_m3s": w["Upstream Flow Arapuni (cumecs, m3/s)"],
        "Q_out_m3s": w["Downstream Flow Karapiro (cumecs, m3/s)"],
    })
    return camp, flows


# ---------------------------------------------------------------------------
# Raw ICP-MS (corrected mode)
# ---------------------------------------------------------------------------
_EXCLUDE_BATCHES = ("Batch_4-_Solution_report_Adam_5.06.25.xlsx",)   # duplicate export


def read_icpms_long(path, exclude_batches=_EXCLUDE_BATCHES):
    """Long ICP-MS table: one row per measured sample.
    Columns needed: site_code, kind (F-CAT|U-MET), date, rep, Ca_ugL, Mg_ugL (batch optional)."""
    d = pd.read_csv(path, parse_dates=["date"])
    if "batch" in d:
        d = d[~d["batch"].isin(exclude_batches)]
    return d


def ca_by_estimator(icp, estimator, site_map=None):
    """Collapse replicate ICP-MS samples to one Ca value (mg/L) per site and date.

    fcat_mean      mean of the filtered (F-CAT) replicates (as published)
    median4        median of all F-CAT and U-MET replicates
    mg_normalised  median Ca/Mg ratio x the date's mean of site-median Mg;
                   removes whole-sample dilution artefacts (a check, not
                   the headline estimator)
    """
    site_map = site_map or _site_code_map()
    d = icp.copy()
    d["site"] = d["site_code"].map(site_map)
    d = d.dropna(subset=["site", "Ca_ugL"])
    if estimator == "fcat_mean":
        g = d[d.kind == "F-CAT"].groupby(["date", "site"]).Ca_ugL.mean()
    elif estimator == "median4":
        g = d.groupby(["date", "site"]).Ca_ugL.median()
    elif estimator == "mg_normalised":
        d["CaMg"] = d.Ca_ugL / d.Mg_ugL
        ratio = d.groupby(["date", "site"]).CaMg.median()
        mg = d.groupby(["date", "site"]).Mg_ugL.median()
        mg_date = mg.groupby(level="date").mean()
        g = ratio * mg_date.reindex(ratio.index.get_level_values("date")).values
    else:
        raise ValueError(f"unknown Ca estimator {estimator!r}")
    return (g / 1e3).rename("Ca_mgL").reset_index()


def _site_code_map():
    # site codes as written on sample labels
    return {"S1": "up", "S3": "down", "S2B": "bot", "S2M": "mid", "S2T": "top",
            "S4": "narrows"}


def read_robust_site_date(path, column="Ca_med4"):
    """Pre-collapsed robust table (data/update/icpms/water_ca_site_date_robust.csv).
    Columns available: Ca_FCAT, Ca_med4, Ca_Mgnorm (ug/L)."""
    d = pd.read_csv(path, parse_dates=["date"])
    d["site"] = d["site"].map(_ROBUST_SITE)
    out = d[["date", "site"]].copy()
    out["Ca_mgL"] = d[column] / 1e3
    return out.dropna()


def read_alkalinity(path):
    """data/update/alk/alkalinity_site_date.csv: date, site, alk45 (mg/L CaCO3)."""
    d = pd.read_csv(path, parse_dates=["date"])
    d["site"] = d["site"].map(_ROBUST_SITE).fillna(d["site"])
    return d.rename(columns={"alk45": "Alk_mgL"})[["date", "site", "Alk_mgL"]]


# ---------------------------------------------------------------------------
# Flows
# ---------------------------------------------------------------------------
def read_daily_flows(path):
    """Daily Karapiro outflow, derived inflow and level.
    Accepts data/update/flow_daily_karapiro_2024_2026.csv (DateTime, Flow, Qin, Level)
    or the template layout (date, Q_out_m3s, Q_in_m3s, level_m)."""
    d = pd.read_csv(path)
    d = d.rename(columns={"DateTime": "date", "Flow": "Q_out_m3s",
                          "Qin": "Q_in_m3s", "Level": "level_m"})
    d["date"] = pd.to_datetime(d["date"])
    return d


def derive_inflow(outflow_15min, level_15min, surface_area_m2):
    """Inflow = outflow + A_surface * dh/dt, daily means.
    Inputs are pandas Series indexed by timestamp (m3/s and m).
    Against the published Arapuni flows this gives r2 = 0.991 at Karapiro."""
    q = outflow_15min.resample("D").mean()
    h = level_15min.resample("D").mean()
    dhdt = h.diff().shift(-1).add(h.diff()).div(2) / 86400     # centred, m/s
    return pd.DataFrame({"Q_out_m3s": q, "Q_in_m3s": q + surface_area_m2 * dhdt,
                         "level_m": h}).rename_axis("date").reset_index()


# ---------------------------------------------------------------------------
# Specimens and surveys
# ---------------------------------------------------------------------------
def read_specimens(path, sheet=None):
    """Specimen table -> sample, dims (mm), dry_total_g, dry_shell_g, volume_cm3 (ellipsoid).

    Reads either the original sheet ('corbicula shell mass', columns 'Dry total',
    'Dry shell') or the corrected sheet from the lab notes ('corbicula shell mass
    corrected', columns dry_total_g, dry_shell_g plus tray and gross weights).
    With no sheet given, the corrected sheet is used when present."""
    xl = pd.ExcelFile(path)
    if sheet is None:
        sheet = ("corbicula shell mass corrected" if "corbicula shell mass corrected" in xl.sheet_names
                 else "corbicula shell mass")
    s = pd.read_excel(xl, sheet_name=sheet)
    s.columns = s.columns.str.strip()
    s = s[s["sample"].notna()].copy()
    total = "dry_total_g" if "dry_total_g" in s else "Dry total"
    shell = "dry_shell_g" if "dry_shell_g" in s else "Dry shell"
    out = pd.DataFrame({
        "sample": s["sample"].astype(str).str.strip(),
        "length_mm": pd.to_numeric(s["length_mm"]),
        "height_mm": pd.to_numeric(s["height_mm"]),
        "width_mm": pd.to_numeric(s["width_mm"]),
        "dry_total_g": pd.to_numeric(s[total]),
        "dry_shell_g": pd.to_numeric(s[shell]),
    })
    if "tray_wt_g" in s and "gross_total_g" in s:
        gap = (pd.to_numeric(s["gross_total_g"]) - pd.to_numeric(s["tray_wt_g"]) - out["dry_total_g"]).abs()
        if (gap > 1e-3).any():
            raise ValueError(f"tare check failed for {list(out.loc[gap > 1e-3, 'sample'])}")
    out["volume_cm3"] = (np.pi / 6) * out.length_mm * out.height_mm * out.width_mm / 1e3
    return out


def read_surveys(path, sheet="corbicula density"):
    s = pd.read_excel(path, sheet_name=sheet)
    s.columns = s.columns.str.strip()
    return s.rename(columns={"mean density_m2": "density_m2"})
