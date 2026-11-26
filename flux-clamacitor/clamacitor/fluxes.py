"""Mass balances: Ca and alkalinity removal as CaCO3 (mol/d and t/d).

All functions take a wide table with one row per sampling date and columns
Ca_<site>, Alk_<site> (mg/L), Q_in_m3s, Q_out_m3s. Positive J = removal.
"""
import numpy as np
import pandas as pd
from .constants import M_CACO3, EQW_ALK, SEC_PER_DAY


def to_wide(campaign, flows):
    """campaign long (date, site, Ca_mgL[, Alk_mgL]) + flows -> wide by date."""
    vals = [c for c in ("Ca_mgL", "Alk_mgL") if c in campaign]
    w = campaign.pivot_table(index="date", columns="site", values=vals)
    w.columns = [f"{v.split('_')[0]}_{s}" for v, s in w.columns]
    f = flows.set_index("date")[["Q_in_m3s", "Q_out_m3s"]]
    return w.join(f, how="left").sort_index()


def longitudinal(w, method, analyte="Ca", lake=None):
    """Eq. 1/2. 'published': (C_up Q_in - C_down Q_out); 'qin_delta': Q_in (C_up - C_down).
    The published form books any inflow/outflow imbalance on the day as removal."""
    up, dn = w[f"{analyte}_up"], w[f"{analyte}_down"]
    if method.longitudinal == "published":
        flux = up * w.Q_in_m3s - dn * w.Q_out_m3s
    elif method.longitudinal == "qin_delta":
        flux = (up - dn) * w.Q_in_m3s
    else:
        raise ValueError(method.longitudinal)
    J = _to_mol(flux, analyte, method)
    if method.storage_term:
        J = J - storage(w, analyte, lake, method)
    return J


def vertical(w, method, analyte="Ca"):
    """Eq. 3. Top-minus-bottom at the dam face x outflow."""
    q = w.Q_out_m3s if method.vertical_flow == "Q_down" else w.Q_in_m3s
    return _to_mol((w[f"{analyte}_top"] - w[f"{analyte}_bot"]) * q, analyte, method)


def storage(w, analyte, lake, method):
    """V dC/dt for the lake (mol CaCO3/d), using the downstream station as the
    lake concentration and centred differences between sampling dates."""
    c = w[f"{analyte}_down"]
    t = (w.index - w.index[0]).days.values.astype(float)
    dcdt = pd.Series(np.gradient(c.values, t), index=w.index)        # mg/L/d
    mass = dcdt * lake.volume_m3                                      # g/d
    if analyte == "Ca":
        return mass / method.m_ca
    return mass / EQW_ALK / 2


def _to_mol(flux_mg_per_L_m3s, analyte, method):
    # mg/L * m3/s = g/s
    g_per_day = flux_mg_per_L_m3s * SEC_PER_DAY
    if analyte == "Ca":
        return g_per_day / method.m_ca
    return g_per_day / EQW_ALK / 2


def t_caco3(J_mol_d):
    return J_mol_d * M_CACO3 * 1e-6


def all_fluxes(w, method, lake):
    """Return a per-date table of removal in t CaCO3/d for the four balances."""
    out = pd.DataFrame(index=w.index)
    out["tau_d"] = lake.volume_m3 / (w.Q_out_m3s * SEC_PER_DAY)
    for a in ("Ca", "Alk"):
        if f"{a}_up" in w and f"{a}_down" in w:
            out[f"t_{a}_L"] = t_caco3(longitudinal(w, method, a, lake))
        if f"{a}_top" in w and f"{a}_bot" in w:
            out[f"t_{a}_V"] = t_caco3(vertical(w, method, a))
    return out


def cumulative(dates, t_per_day, clip_negative=False):
    """Trapezoidal integral of t/d between sampling dates (t CaCO3)."""
    s = pd.Series(np.asarray(t_per_day, float), index=pd.DatetimeIndex(dates)).dropna()
    y = s.clip(lower=0).values if clip_negative else s.values
    days = (s.index - s.index[0]).days.values
    return pd.Series(np.concatenate([[0.0], np.cumsum(np.diff(days) * (y[1:] + y[:-1]) / 2)]),
                     index=s.index)
