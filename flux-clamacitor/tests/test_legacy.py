"""Legacy mode must reproduce Hartland et al. (2026) Table 1 / Section 3.4."""
import pytest
from conftest import PUBLISHED, need
import clamacitor as cc
from clamacitor import scaling

pytestmark = need(PUBLISHED)


@pytest.fixture(scope="module")
def r():
    return cc.run_legacy(PUBLISHED)


def test_sampling_dates(r):
    assert r.summary["n_dates"] == 17
    assert r.summary["first_date"] == "2024-10-17"
    assert r.summary["last_date"] == "2025-07-23"


def test_ca_vertical_biovolume(r):
    b = r.summary["balances"]["Ca_V"]
    assert b["n_pos"] == 13
    assert round(b["Vbio_mean"], 1) == 14.0
    assert round(b["Vbio_median"], 1) == 8.7


def test_removal_rate_abstract(r):
    # abstract: ~14 +/- 14 t/d
    b = r.summary["balances"]["Ca_V"]
    assert round(b["t_d_mean"]) == 14 and round(b["t_d_sd"]) == 14


@pytest.mark.parametrize("bal,mean,sd,N", [("Ca_V", 343, 331, 2.54), ("Alk_V", 401, 342, 2.97)])
def test_density_by_method(r, bal, mean, sd, N):
    b = r.summary["balances"][bal]
    assert round(b["D_mean"]) == mean
    assert round(b["D_sd"]) == sd
    assert round(b["N_mean"] / 1e9, 2) == N


def test_combined_density(r):
    p = r.summary["pooled"]
    assert p["n"] == 24
    assert round(p["D_mean"]) == 370
    assert round(p["D_sd"]) == 323
    assert round(p["D_median"]) == 358
    assert round(p["N_mean"] / 1e9, 2) == 2.74
    assert round(p["Vbio_mean"], 1) == 15.1


def test_published_occupancy_arithmetic(r):
    # N / survey mean density. The paper states 3.94e6 m2 (63 %); the
    # arithmetic gives 4.92e6 m2 (66 %). Locked in here as a correction item.
    area = scaling.occupied_area(r.summary["pooled"]["N_mean"], 557)
    assert round(float(area) / 1e6, 2) == 4.92
