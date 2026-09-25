#!/usr/bin/env python
"""R1.4: minimum detectable effect (MDE) of each fixed-panel randomization test.

Scope: these MDE values are marginal normal-shift sensitivity approximations,
conditional on the frozen panel and the named null reference distributions. They
are not guaranteed joint 80%-power thresholds, population power calculations, or
exclusion bounds on an unobserved effect or regulatory capability.

Two versions:
  analytic   MDE = (z_{1-a/2} + z_{0.80}) * sd_null
             assumes a centered normal-shift reference
  empirical  MDE = c_a + z_{0.80} * sd_null, c_a the (1-a) quantile of |null|
             from saved N = 9,999 arrays; still an approximate marginal shift model
with a = 0.05 per test and a = 0.05/m (a fixed Bonferroni reference, not adaptive BH
family power). Historical fields named mde_dual contain the maximum of two
marginal approximations; taking this maximum does not ensure joint 80% power.
power_sim.py evaluates joint rejection only for its specified synthetic mechanism;
it does not validate the assumptions for every model-specific null distribution.
"""
import argparse
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary-only", action="store_true",
                    help="Use only published null standard deviations; no empirical null arrays")
    ap.add_argument("--output", type=Path,
                    default=Path(common.fpa.ROOT) / "output/revision_checks/mde.json")
    args = ap.parse_args()
    doc = json.loads(common.AUDIT_JSON.read_text())
    res_dir = common.OUT_DIR / "resolution_n9999"
    have_deep = not args.summary_only and res_dir.exists() and any(res_dir.glob("fm_*.npz"))
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
    probe_path = common.result_json("tf_probe_pair_stats_v2.json", "SCREG_PROBE_JSON")
    probe = json.loads(probe_path.read_text())
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
        "scope": ("Marginal normal-shift sensitivity approximations on a frozen panel; "
                  "mde_dual fields are maximum-marginal references, not guaranteed joint "
                  "80%-power thresholds, population power, or effect-size exclusion bounds."),
        "power": POWER, "alpha": ALPHA,
        "rows": rows_out,
        "summary_ranges": summary,
        "probe": probe_out,
        "inputs": {"fixed_panel_audit_v2_sha256": common.sha256_file(common.AUDIT_JSON),
                   "tf_probe_pair_stats_v2_sha256": common.sha256_file(probe_path),
                   "n9999_nulls_used": bool(deep)},
    }
    common.write_json(args.output, result)
    for k, v in summary.items():
        log(k, {kk: [round(x, 4) for x in vv] for kk, vv in v.items()})


if __name__ == "__main__":
    main()
