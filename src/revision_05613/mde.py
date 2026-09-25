#!/usr/bin/env python
"""R1.4: minimum detectable effect (MDE) of each fixed-panel randomization test.

Scope: the MDE is the smallest partial Spearman rho that the named randomization
test would declare unusual with 80% probability on this frozen panel. It is a
statement about the test's sensitivity under the named randomizations, not a
population power calculation and not an exclusion bound on regulatory capability.

Two versions:
  analytic   MDE = (z_{1-a/2} + z_{0.80}) * sd_null           (normal approximation)
  empirical  MDE = c_a + z_{0.80} * sd_null, c_a the (1-a) quantile of |null|
             from the saved N = 9,999 null arrays (when present)
with a = 0.05 per test and a = 0.05/m (Bonferroni bound on BH's first rejection).
The dual-null MDE of a row is the larger of its two single-null MDEs.
power_sim.py checks the normal approximation by simulation.
"""
import json
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import log  # noqa: E402

Z = NormalDist().inv_cdf
POWER = 0.80
ALPHA = 0.05


def analytic(sd, a):
    return (Z(1 - a / 2) + Z(POWER)) * sd


def main():
    doc = json.loads(common.AUDIT_JSON.read_text())
    res_dir = common.OUT_DIR / "resolution_n9999"
    have_deep = res_dir.exists() and any(res_dir.glob("fm_*.npz"))
    deep = {}
    if have_deep:
        # file names: fm_{group}_{null}_{spec}.npz, where group and spec may contain '_'
        for f in res_dir.glob("fm_*.npz"):
            stem = f.stem[3:]
            spec = next((s for s in ("non_degree", "full") if stem.endswith("_" + s)), None)
            if spec is None:
                continue
            group, null_type = stem[: -len(spec) - 1].rsplit("_", 1)
            Zf = np.load(f)
            for lab in Zf.files:
                if not lab.startswith("_"):
                    tissue = "brain" if group == "brain" else "pbmc"
                    deep[(tissue, spec, lab, null_type)] = Zf[lab]
    rows_out = []
    for tissue in ("brain", "pbmc"):
        for fam in ("primary_family", "sensitivity_family"):
            rows = doc["pooled"][tissue][fam]["rows"]
            m = len(rows)
            for r in rows:
                spec = r["confound_spec"]
                rec = {"tissue": tissue, "spec": spec, "model_label": r["model_label"],
                       "rho": r["observed_partial_rho"], "m_family": m}
                for nt, key in (("mantel", "mantel"), ("degree", "degree_preserving")):
                    sd = r[key]["null_sd"]
                    rec[f"sd_{nt}_n999"] = sd
                    rec[f"mde_{nt}_a05"] = analytic(sd, ALPHA)
                    rec[f"mde_{nt}_a05m"] = analytic(sd, ALPHA / m)
                    arr = deep.get((tissue, spec, r["model_label"], nt))
                    if arr is not None:
                        a = np.abs(arr)
                        rec[f"sd_{nt}_n9999"] = float(arr.std())
                        for tag, lev in (("a05", ALPHA), ("a05m", ALPHA / m)):
                            c = float(np.quantile(a, 1 - lev))
                            rec[f"crit_{nt}_{tag}_emp"] = c
                            rec[f"mde_{nt}_{tag}_emp"] = c + Z(POWER) * float(arr.std())
                for tag in ("a05", "a05m"):
                    rec[f"mde_dual_{tag}"] = max(rec[f"mde_mantel_{tag}"], rec[f"mde_degree_{tag}"])
                    if f"mde_mantel_{tag}_emp" in rec:
                        rec[f"mde_dual_{tag}_emp"] = max(rec[f"mde_mantel_{tag}_emp"],
                                                         rec[f"mde_degree_{tag}_emp"])
                rec["abs_rho_below_dual_mde_a05m"] = abs(rec["rho"]) < rec["mde_dual_a05m"]
                rows_out.append(rec)

    # TF-disjoint probe (Table 8): family permutation and paired sign-flip contrasts.
    probe = json.loads((Path(common.fpa.OUT) / "tf_probe_pair_stats_v2.json").read_text())
    probe_out = {"families": {}, "contrasts": {}}
    for fam, d in probe["families"].items():
        probe_out["families"][fam] = {"rho": d["adjusted_rho_mean"], "sd": d["null_std"],
                                      "mde_a05": analytic(d["null_std"], ALPHA)}
    for fam, d in probe["contrasts_vs_baseline"].items():
        se = d["paired_delta_std"] / np.sqrt(d["n_tf"])
        probe_out["contrasts"][fam] = {"delta": d["paired_delta_mean"], "se": float(se),
                                       "mde_a05": analytic(float(se), ALPHA)}

    def rng_(sel, key):
        v = [r[key] for r in rows_out if sel(r)]
        return [float(min(v)), float(max(v))] if v else None

    summary = {}
    for tissue in ("brain", "pbmc"):
        for spec in ("full", "non_degree"):
            sel = lambda r, t=tissue, s=spec: r["tissue"] == t and r["spec"] == s  # noqa: E731
            summary[f"{tissue}_{spec}"] = {k: rng_(sel, k) for k in (
                "sd_mantel_n999", "sd_degree_n999", "mde_dual_a05", "mde_dual_a05m")}
    result = {
        "schema_version": 1,
        "analysis": "minimum_detectable_effect",
        "scope": ("Sensitivity of the fixed-panel randomization tests under the named "
                  "randomizations; not population power, not an exclusion bound on "
                  "regulatory capability."),
        "power": POWER, "alpha": ALPHA,
        "rows": rows_out,
        "summary_ranges": summary,
        "probe": probe_out,
        "inputs": {"fixed_panel_audit_v2_sha256": common.sha256_file(common.AUDIT_JSON),
                   "n9999_nulls_used": bool(deep)},
    }
    common.write_json(common.OUT_DIR / "mde.json", result)
    for k, v in summary.items():
        log(k, {kk: [round(x, 4) for x in vv] for kk, vv in v.items()})


if __name__ == "__main__":
    main()
