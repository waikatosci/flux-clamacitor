import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
DATA = ROOT / "data"
PUBLISHED = DATA / "Full_data_for_Corbicula_Calcs_Landscape_Updated_16_1_26.xlsx"
ICPMS = DATA / "correction" / "icpms_replicates_2024-10_2025-07.csv"
SPEC = DATA / "Density and shell mass for AH.xlsx"                          # as published
SPEC2 = DATA / "correction" / "shell_mass_corrected_2026-11-26.xlsx"        # corrected weighings


def need(*paths):
    missing = [p.name for p in paths if not p.exists()]
    return pytest.mark.skipif(bool(missing), reason=f"missing data: {missing}")
