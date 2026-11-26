"""Flux Clamacitor 2.0: CaCO3 removal mass balances and population inversion
for invasive bivalves (Hartland et al. 2026, ES&T, doi:10.1021/acs.est.6c03421)."""
from .constants import Lake, Method, LEGACY, CORRECTION, PRESETS
from .pipeline import run_legacy, run_correction, run_tables

__version__ = "2.0.0"
