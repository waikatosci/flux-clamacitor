"""Command line:  python -m clamacitor {legacy|correction|tables} ...

  legacy      reproduce the published ES&T values from the published workbook
  correction  the values in the Addition and Correction
  tables      run on the input templates (campaign + flows) for any lake
"""
import argparse, json, sys
from pathlib import Path
from . import config, pipeline
from .constants import Lake, PRESETS


def _round(o):
    if isinstance(o, float):
        return float(f"{o:.4g}")
    if isinstance(o, dict):
        return {k: _round(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_round(v) for v in o]
    return o


def main(argv=None):
    ap = argparse.ArgumentParser(prog="clamacitor", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["legacy", "correction", "tables"])
    ap.add_argument("--config", default="configs/karapiro.yaml")
    ap.add_argument("--campaign", help="campaign template CSV (tables mode)")
    ap.add_argument("--flows", help="flows CSV (overrides config)")
    ap.add_argument("--preset", choices=list(PRESETS), help="method preset (tables mode)")
    ap.add_argument("--start"); ap.add_argument("--end")
    ap.add_argument("--out", help="folder for per_date.csv and summary.json")
    a = ap.parse_args(argv)

    lake, method, inp, period = config.load(a.config)
    start = a.start or period.get("start"); end = a.end or period.get("end")
    if a.mode == "legacy":
        r = pipeline.run_legacy(inp["published"], lake=lake)
    elif a.mode == "correction":
        r = pipeline.run_correction(inp["published"], inp["icpms"], lake=lake, method=method)
    else:
        if not a.campaign:
            ap.error("tables mode needs --campaign")
        m = PRESETS[a.preset] if a.preset else method
        if not a.flows:
            ap.error("tables mode needs --flows")
        r = pipeline.run_tables(a.campaign, a.flows, lake=lake, method=m,
                                start=start, end=end)
    s = _round(r.summary)
    print(json.dumps(s, indent=2))
    if a.out:
        o = Path(a.out); o.mkdir(parents=True, exist_ok=True)
        r.per_date.to_csv(o / "per_date.csv")
        (o / "summary.json").write_text(json.dumps(s, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
