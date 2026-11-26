"""End-to-end runs for the two presets, returning per-date tables and a summary."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd

from . import io, fluxes, scaling
from .constants import Lake, Method, LEGACY, CORRECTION

BALANCES = ("Ca_V", "Alk_V", "Ca_L", "Alk_L")


@dataclass
class Result:
    method: Method
    lake: Lake
    per_date: pd.DataFrame       # t CaCO3/d, biovolume, N, D per balance
    summary: dict


# ---------------------------------------------------------------------------
def run_legacy(published_xlsx, lake=Lake(), method=LEGACY):
    camp, flows = io.read_published_wide(published_xlsx)
    w = fluxes.to_wide(camp, flows)
    return _finish(w, method, lake, balances=("Ca_V", "Alk_V", "Ca_L", "Alk_L"),
                   pooled=("Ca_V", "Alk_V"))


def run_correction(published_xlsx, icpms_csv, lake=Lake(), method=CORRECTION):
    """Values in the Addition and Correction: the published workbook, flows and
    equations, with Ca replaced by the robust replicate estimate."""
    camp, flows = io.read_published_wide(published_xlsx)
    ca = io.ca_by_estimator(io.read_icpms_long(icpms_csv), method.ca_estimator)
    camp = camp.drop(columns="Ca_mgL").merge(ca, on=["date", "site"], how="left")
    w = fluxes.to_wide(camp, flows)
    return _finish(w, method, lake, balances=("Ca_V", "Alk_V", "Ca_L", "Alk_L"),
                   pooled=("Ca_V", "Alk_V"))


# ---------------------------------------------------------------------------
def _finish(w, method, lake, balances, pooled):
    f = fluxes.all_fluxes(w, method, lake)
    for b in balances:
        if f"t_{b}" not in f:
            continue
        f[f"Vbio_{b}"] = scaling.biovolume(f[f"t_{b}"], method)
        f[f"N_{b}"] = scaling.population(f[f"Vbio_{b}"], method)
        f[f"D_{b}"] = scaling.density(f[f"N_{b}"], lake)
    summ = {"method": method.name, "n_dates": int(len(f)),
            "first_date": str(f.index.min().date()), "last_date": str(f.index.max().date()),
            "days": int((f.index.max() - f.index.min()).days),
            "balances": {}}
    for b in balances:
        if f"t_{b}" not in f:
            continue
        summ["balances"][b] = _stats(f, b, method, lake)
    summ["pooled"] = _pooled(f, pooled, method)
    return Result(method, lake, f, summ)


def _stats(f, b, method, lake):
    t = f[f"t_{b}"].dropna()
    pos = t[t > 0]
    use = pos if method.averaging == "positive_only" else t
    out = {"n": int(len(t)), "n_pos": int(len(pos)),
           "t_d_mean": float(use.mean()), "t_d_sd": _sd(use),
           "t_d_median": float(use.median())}
    out.update(_scale_block(use, method, lake))
    out["cum_t"] = float(fluxes.cumulative(t.index, t).iloc[-1])
    if method.report_bound and len(pos):
        out["upper_bound_positive_only"] = {"t_d_mean": float(pos.mean()),
                                            **_scale_block(pos, method, lake)}
    if method.bulk_density_range:
        lo, hi = method.bulk_density_range
        out["D_mean_range"] = [float(_scale_block(use, method, lake, rho=hi)["D_mean"]),
                               float(_scale_block(use, method, lake, rho=lo)["D_mean"])]
    return out


def _sd(x, ddof=1):
    x = np.asarray(x, float)
    return float(np.std(x, ddof=ddof)) if x.size > ddof else None


def _scale_block(t, method, lake, rho=None):
    v = scaling.biovolume(t, method, rho)
    N = scaling.population(v, method)
    D = scaling.density(N, lake)
    return {"Vbio_mean": float(np.mean(v)), "Vbio_median": float(np.median(v)),
            "N_mean": float(np.mean(N)), "D_mean": float(np.mean(D)),
            "D_sd": _sd(D), "D_median": float(np.median(D))}


def _pooled(f, pooled, method):
    """Pool dates across balances, as in the published combined figure."""
    parts = []
    for b in pooled:
        if f"D_{b}" in f:
            d = f[[f"t_{b}", f"D_{b}", f"Vbio_{b}", f"N_{b}"]].dropna()
            d.columns = ["t", "D", "V", "N"]
            if method.averaging == "positive_only":
                d = d[d.t > 0]
            parts.append(d)
    if not parts:
        return {}
    p = pd.concat(parts)
    return {"balances": list(pooled), "n": int(len(p)),
            "D_mean": float(p.D.mean()), "D_sd": _sd(p.D, method.sd_ddof_pooled),
            "D_median": float(p.D.median()), "Vbio_mean": float(p.V.mean()),
            "N_mean": float(p.N.mean())}


# ---------------------------------------------------------------------------
def run_tables(campaign_csv, flows_csv, lake=Lake(), method=CORRECTION, start=None, end=None):
    """Generic entry point for the input templates (and the web estimator).

    campaign_csv: date, site, Ca_mgL[, Alk_mgL]  (one value per site and date,
                  already collapsed from replicates). Sites: up, down, top, bot.
    flows_csv:    date, Q_out_m3s, Q_in_m3s
    """
    camp = pd.read_csv(campaign_csv, parse_dates=["date"])
    flows = io.read_daily_flows(flows_csv)
    w = fluxes.to_wide(camp, flows)
    if start:
        w = w[w.index >= pd.Timestamp(start)]
    if end:
        w = w[w.index <= pd.Timestamp(end)]
    w = w.dropna(subset=["Q_in_m3s", "Q_out_m3s"])
    bal = tuple(b for b in BALANCES
                if {f"{b.split('_')[0]}_{s}" for s in (("top", "bot") if b.endswith("V") else ("up", "down"))} <= set(w.columns))
    return _finish(w, method, lake, balances=bal, pooled=("Ca_V",))
