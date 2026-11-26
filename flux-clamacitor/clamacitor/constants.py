"""Physical constants and the two method presets.

`LEGACY` reproduces Hartland et al. (2026, ES&T, doi:10.1021/acs.est.6c03421)
exactly. `CORRECTION` gives the values in the Addition and Correction to that
article. Any field can be overridden from a lake config (YAML).
"""
from dataclasses import dataclass, field, replace, asdict

M_CA = 40.078          # g/mol (published scripts used 40.08; difference < 0.01 %)
M_CA_LEGACY = 40.08
M_CACO3 = 100.09       # g/mol
EQW_ALK = 50.043       # g CaCO3 per eq (alkalinity reported as mg/L CaCO3)
SEC_PER_DAY = 86400.0


@dataclass(frozen=True)
class Lake:
    name: str = "Karapiro"
    volume_m3: float = 85e6
    benthic_area_m2: float = 7.409e6      # area used to scale N to density
    surface_area_m2: float = 7.4e6        # used to derive inflow from level change


@dataclass(frozen=True)
class Method:
    name: str
    # chemistry
    ca_estimator: str            # 'fcat_mean' | 'median4' | 'mg_normalised'
    # flux equation for the longitudinal balance
    longitudinal: str            # 'published' (CinQin - CoutQout) | 'qin_delta'
    storage_term: bool           # subtract V dC/dt (longitudinal only)
    vertical_flow: str           # 'Q_down' (outflow) – same in both presets
    # averaging across sampling dates
    averaging: str               # 'positive_only' | 'all_dates'
    report_bound: bool           # also report the positive-only figure as an upper bound
    # biometrics
    total_shell_ratio: float = 1.13       # g total dry mass / g shell CaCO3
    bulk_density_t_m3: float = 1.15       # g/cm3 == t/m3
    bulk_density_range: tuple = ()        # (lo, hi) for sensitivity, optional
    dvdt_cm3_d: float = 0.0055            # per-individual biovolume growth
    sd_ddof_pooled: int = 0               # published pooled SD used np.std (ddof=0)
    m_ca: float = M_CA_LEGACY

    def with_(self, **kw):
        return replace(self, **kw)

    def as_dict(self):
        return asdict(self)


LEGACY = Method(
    name="legacy",
    ca_estimator="fcat_mean",
    longitudinal="published",
    storage_term=False,
    vertical_flow="Q_down",
    averaging="positive_only",
    report_bound=False,
)

CORRECTION = Method(
    # Addition and Correction to Hartland et al. (2026): the published method,
    # with replicate Ca taken as the median of the four filtered and unfiltered
    # measurements and the bulk density from the corrected specimen weighings.
    name="correction",
    ca_estimator="median4",
    longitudinal="published",
    storage_term=False,
    vertical_flow="Q_down",
    averaging="positive_only",
    report_bound=False,
    bulk_density_t_m3=0.632,             # total dry mass / total volume, 29 specimens
    bulk_density_range=(0.57, 0.70),     # bootstrap 95% interval over specimens
)

PRESETS = {"legacy": LEGACY, "correction": CORRECTION}
