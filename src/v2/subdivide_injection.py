#!/usr/bin/env python
"""Subdivide the injection ladder at alpha in {0.002, 0.005, 0.01}.

effect_vs_injection_scale_v2.json documents that every full-spec pooled effect in
the paper is smaller than the rho recovered at alpha=0.02, the smallest non-zero
injection probed -- the ladder has no calibration point at or below the reported
effects' magnitude. This runs the audit's own run_sensitivity with three finer
alphas (30 replicates each, both tissues) and places the observed pooled effects
on the subdivided curve.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, ".")
import fixed_panel_audit as fpa  # noqa: E402
import run_fixed_panel_audit as R  # noqa: E402

OUT = fpa.OUT
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
ATAC_B = os.environ.get(
    "SCFM_BRAIN_ATAC",
    f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")
ATAC_P = f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"
FINE_ALPHAS = [0.002, 0.005, 0.01]
N_REP = 30


def main():
    R.ALPHAS = FINE_ALPHAS  # run_sensitivity reads the module constant

    G_b, co_b, _, tf_b, _ = R.load_pooled_brain()
    G_p, co_p, _, tf_p, _ = R.load_pooled_pbmc()

    ss = np.random.SeedSequence([R.SEED_ROOT, 20260730, 7])
    out = {}
    for tissue, G, co, tf, atac in (
        ("brain", G_b, co_b, tf_b, ATAC_B),
        ("pbmc", G_p, co_p, tf_p, ATAC_P),
    ):
        print(f"--- {tissue} ---", flush=True)
        res = R.run_sensitivity(ss.spawn(1)[0], tissue, G, co, tf, atac, N_REP)
        out[tissue] = res
        for row in res["rows"]:
            vals = [r["observed_partial_rho_axis_aligned"] for r in row["replicate_runs"]]
            print(f"  alpha={row['alpha']:.3f}  mean rho={np.mean(vals):.6f} "
                  f"range=[{min(vals):.6f},{max(vals):.6f}]", flush=True)

    # Place the observed pooled full-spec effects on the curve
    audit = json.load(open(f"{OUT}/fixed_panel_audit_v2.json"))
    effects = []
    for tissue in ("brain", "pbmc"):
        for r in audit["pooled"][tissue]["primary_family"]["rows"]:
            effects.append({
                "tissue": tissue, "model_label": r["model_label"],
                "observed_rho": r["observed_partial_rho"],
            })
    out["observed_effects"] = effects
    with open(f"{OUT}/injection_subdivided_v2.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote injection_subdivided_v2.json", flush=True)


if __name__ == "__main__":
    main()
