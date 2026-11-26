"""Template pathway: the same code the web estimator will call."""
import pandas as pd
import clamacitor as cc
from clamacitor import pipeline


def test_tables_roundtrip(tmp_path):
    camp = pd.DataFrame({
        "date": ["2025-01-10"] * 4 + ["2025-01-24"] * 4,
        "site": ["up", "down", "top", "bot"] * 2,
        "Ca_mgL": [6.0, 5.6, 5.4, 5.2, 6.1, 5.8, 5.5, 5.4],
    })
    flows = pd.DataFrame({"date": ["2025-01-10", "2025-01-24"],
                          "Q_out_m3s": [150.0, 150.0], "Q_in_m3s": [150.0, 150.0]})
    camp.to_csv(tmp_path / "c.csv", index=False)
    flows.to_csv(tmp_path / "f.csv", index=False)
    r = pipeline.run_tables(tmp_path / "c.csv", tmp_path / "f.csv")
    # 0.2 and 0.1 mg/L x 150 m3/s -> mol Ca/d -> t CaCO3/d
    exp = [0.2 * 150 * 86400 / 40.08 * 100.09e-6, 0.1 * 150 * 86400 / 40.08 * 100.09e-6]
    assert list(r.per_date["t_Ca_V"].round(6)) == [round(x, 6) for x in exp]
    assert set(r.summary["balances"]) == {"Ca_V", "Ca_L"}
