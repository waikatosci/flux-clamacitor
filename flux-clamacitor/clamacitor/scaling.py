"""Removal -> biovolume -> population -> density, and the occupied-area check."""
import numpy as np
import pandas as pd


def biovolume(t_per_day, method, bulk_density=None):
    """Eq. 5. m3 of live biomass produced per day."""
    rho = method.bulk_density_t_m3 if bulk_density is None else bulk_density
    return np.asarray(t_per_day, float) * method.total_shell_ratio / rho


def population(vbio_m3_d, method):
    """Eq. 6. N = V_bio / (dV/dt)."""
    return np.asarray(vbio_m3_d, float) / (method.dvdt_cm3_d * 1e-6)


def density(N, lake):
    """Eq. 7. individuals per m2 over the benthic area."""
    return np.asarray(N, float) / lake.benthic_area_m2


def occupied_area(N, survey_density):
    """Area (m2) the population would occupy at the surveyed density."""
    return np.asarray(N, float) / survey_density


def bulk_density_from_specimens(spec, exclude=()):
    """Total dry mass / total ellipsoid volume (g/cm3).

    A ratio of totals weights each specimen by its size. The published 1.15
    was a mean of per-specimen ratios, taken before transcription and tare
    errors in three specimens (small1, small3, medium2) were corrected from
    the lab notes (Nov 2026). With the corrected weighings, all 29: 0.632."""
    s = spec[~spec["sample"].isin(exclude)]
    return float(s.dry_total_g.sum() / s.volume_cm3.sum())


def total_shell_ratio_from_specimens(spec, exclude=()):
    s = spec[~spec["sample"].isin(exclude)]
    return float((s.dry_total_g / s.dry_shell_g).mean())


def bulk_density_bootstrap(spec, n=20000, seed=1, q=(2.5, 97.5)):
    """Percentile interval of the ratio of totals, resampling specimens."""
    rng = np.random.default_rng(seed)
    m, v = spec.dry_total_g.values, spec.volume_cm3.values
    idx = rng.integers(0, len(m), (n, len(m)))
    r = m[idx].sum(1) / v[idx].sum(1)
    return tuple(float(x) for x in np.percentile(r, q))


def dvdt_from_specimens(spec, days=300):
    """Median ellipsoid volume / time to reach it (cm3/ind/d). Published: 1.659/300."""
    return float(spec.volume_cm3.median() / days)


def monte_carlo(t_dates, method, lake, n=20000, rho_range=None,
                area_rel_sd=0.0, seed=42):
    """t_dates: removal (t CaCO3/d) on each sampling date.
    Bootstrap sampling dates, draw bulk density uniformly across rho_range
    and (optionally) perturb area. Returns DataFrame of N and D draws."""
    rng = np.random.default_rng(seed)
    v = np.asarray(t_dates, float)
    v = v[np.isfinite(v)]
    t = rng.choice(v, n)                                    # t CaCO3/d per draw
    lo, hi = rho_range or method.bulk_density_range or (method.bulk_density_t_m3,) * 2
    rho = rng.uniform(lo, hi, n)
    vb = t * method.total_shell_ratio / rho
    N = vb / (method.dvdt_cm3_d * 1e-6)
    area = lake.benthic_area_m2 * (1 + area_rel_sd * rng.standard_normal(n))
    return pd.DataFrame({"t_d": t, "rho": rho, "N": N, "D": N / area})
