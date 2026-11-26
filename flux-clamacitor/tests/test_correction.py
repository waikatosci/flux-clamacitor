"""The correction preset must reproduce the values in the Addition and Correction."""
import pytest
from conftest import PUBLISHED, ICPMS, SPEC, SPEC2, need
import clamacitor as cc
from clamacitor import io, scaling, fluxes


@pytest.fixture(scope="module")
def r():
    return cc.run_correction(PUBLISHED, ICPMS)


@need(PUBLISHED, ICPMS)
def test_removal_rate(r):
    b = r.summary["balances"]["Ca_V"]
    assert b["n_pos"] == 13
    assert round(b["t_d_mean"], 1) == 9.1 and round(b["t_d_sd"], 1) == 7.8


@need(PUBLISHED, ICPMS)
@pytest.mark.parametrize("bal,mean,sd,N", [("Ca_V", 398, 340, 2.95), ("Alk_V", 730, 622, 5.41)])
def test_table1_by_basis(r, bal, mean, sd, N):
    b = r.summary["balances"][bal]
    assert round(b["D_mean"]) == mean and round(b["D_sd"]) == sd
    assert round(b["N_mean"] / 1e9, 2) == N


@need(PUBLISHED, ICPMS)
def test_table1_combined(r):
    p = r.summary["pooled"]
    assert p["n"] == 24
    assert round(p["D_mean"]) == 550 and round(p["D_sd"]) == 496
    assert round(p["D_median"]) == 542
    assert round(p["N_mean"] / 1e9, 2) == 4.07
    assert round(p["Vbio_mean"], 1) == 22.4


@need(PUBLISHED, ICPMS)
def test_occupied_area(r):
    area = float(scaling.occupied_area(r.summary["pooled"]["N_mean"], 557))
    assert round(area / 1e6, 1) == 7.3
    assert round(100 * area / 7.409e6) == 99


@need(PUBLISHED, ICPMS)
def test_residence_time_correlation(r):
    from scipy import stats
    f = r.per_date[r.per_date.t_Ca_V > 0]
    rr, p = stats.pearsonr(f.tau_d, f.t_Ca_V)
    assert round(rr, 2) == -0.71 and round(p, 3) == 0.006


@need(PUBLISHED, ICPMS)
def test_replicate_fix_halves_scatter():
    icp = io.read_icpms_long(ICPMS)
    for est, sd in (("fcat_mean", 0.40), ("median4", 0.19)):
        w = io.ca_by_estimator(icp, est).pivot(index="date", columns="site", values="Ca_mgL")
        assert round(float((w.top - w.bot).std()), 2) == sd


@need(SPEC2)
def test_bulk_density_corrected_weighings():
    s = io.read_specimens(SPEC2)          # corrected sheet; tare check passes
    assert len(s) == 29
    assert round(scaling.bulk_density_from_specimens(s), 3) == 0.632
    lo, hi = scaling.bulk_density_bootstrap(s)
    assert round(lo, 2) == 0.57 and round(hi, 2) == 0.70
    assert round(float((s.dry_total_g / s.volume_cm3).mean()), 2) == 0.87
    assert round(scaling.dvdt_from_specimens(s), 4) == 0.0055


@need(SPEC)
def test_published_bulk_density():
    s = io.read_specimens(SPEC, sheet="corbicula shell mass")
    assert round(float((s.dry_total_g / s.volume_cm3).mean()), 2) == 1.15


def test_published_equation_books_flow_imbalance():
    """With identical Ca up and down, the published longitudinal equation still
    returns removal whenever inflow exceeds outflow."""
    import pandas as pd
    w = pd.DataFrame({"Ca_up": [5.0], "Ca_down": [5.0], "Q_in_m3s": [137.0], "Q_out_m3s": [130.0]},
                     index=pd.to_datetime(["2025-01-01"]))
    assert fluxes.t_caco3(fluxes.longitudinal(w, cc.LEGACY)).iloc[0] == pytest.approx(7.5, abs=0.1)
