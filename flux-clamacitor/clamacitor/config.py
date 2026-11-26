"""Lake and method configuration from YAML."""
from pathlib import Path
import yaml
from .constants import Lake, PRESETS


def load(path):
    cfg = yaml.safe_load(Path(path).read_text())
    lk = dict(cfg.get("lake", {}))
    for k, v in lk.items():                 # YAML 1.1 reads '7.4e6' as a string
        if k != "name":
            lk[k] = float(v)
    lake = Lake(**lk)
    method = PRESETS[cfg.get("preset", "correction")]
    over = cfg.get("method", {}) or {}
    for k, v in list(over.items()):
        if isinstance(v, str):
            try:
                over[k] = float(v)
            except ValueError:
                pass
    if "bulk_density_range" in over:
        over["bulk_density_range"] = tuple(float(x) for x in over["bulk_density_range"])
    period = {k: str(v) for k, v in (cfg.get("period") or {}).items()}
    return lake, method.with_(**over), cfg.get("inputs", {}), period
