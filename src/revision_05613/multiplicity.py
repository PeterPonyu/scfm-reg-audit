#!/usr/bin/env python
"""R1.2: are the dual-null Support decisions robust to the dependence between rows?

Inputs: frozen N = 999 p-values in the configured/local audit JSON, falling back
to results/fixed_panel_audit_v2.public.json on a fresh checkout; and,
when present, the N = 9,999 null arrays from resolution_check.py.

Procedures, each applied within the frozen families (tissue x specification x null):
  BH    Benjamini-Hochberg (the registered rule; valid under independence or PRDS)
  BY    Benjamini-Yekutieli (FDR under arbitrary dependence)
  Holm  Holm step-down (FWER under arbitrary dependence)
  IUT   per-row intersection-union p = max(pM, pD), then BH or BY within tissue x spec
  WY    Westfall-Young step-down minP on the joint randomization distribution
        (needs the saved null arrays; rows that share a randomization batch are
        paired replicate by replicate, independent batches are paired by index)

Corrections on each component family followed by intersection are descriptive
screens; they do not automatically control the rowwise conjunction error rate.
IUT-BY supplies within-family conjunction FDR control under arbitrary dependence,
conditional on valid component p-values. The joint WY construction is exploratory:
subset pivotality and the cross-stream pairing assumptions are not established.

BH is cross-checked with SciPy; BY uses SciPy and Holm uses NumPy.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import log  # noqa: E402

ALPHA = 0.05
NULL_KEYS = {"mantel": "mantel", "degree": "degree_preserving"}


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = p[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(q, 1.0)
    return out


def adjust(p, method):
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be a finite one-dimensional array in [0, 1]")
    if method not in {"BH", "BY", "Holm"}:
        raise ValueError(f"Unknown correction: {method}")
    if not len(p):
        return p.copy()
    if method == "BH":
        ours = bh(p)
        ref = false_discovery_control(p, method="bh")
        if not np.allclose(ours, ref, atol=1e-12):
            raise ValueError("Independent BH implementations disagree")
        return ours
    if method == "BY":
        return false_discovery_control(p, method="by")
    order = np.argsort(p, kind="stable")
    corrected = np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1))
    out = np.empty_like(p)
    out[order] = np.minimum(corrected, 1.0)
    return out


def wy_stepdown_minp(null_mat, observed):
    """Westfall-Young step-down minP. null_mat: (N, m) null statistics; observed: (m,)."""
    N, m = null_mat.shape
    a = np.abs(null_mat)
    obs = np.abs(np.asarray(observed))
    raw = (np.sum(a >= obs, axis=0) + 1) / (N + 1)
    # null p-value of every replicate within its own row's null distribution
    pstar = np.empty_like(a)
    for r in range(m):
        srt = np.sort(a[:, r])
        pstar[:, r] = (N - np.searchsorted(srt, a[:, r], side="left")) / N
    order = np.argsort(raw, kind="stable")
    succ_min = np.minimum.accumulate(pstar[:, order][:, ::-1], axis=1)[:, ::-1]
    adj_sorted = (np.sum(succ_min <= raw[order][None, :], axis=0) + 1) / (N + 1)
    adj_sorted = np.maximum.accumulate(adj_sorted)
    adj = np.empty(m)
    adj[order] = np.minimum(adj_sorted, 1.0)
    return raw, adj


def family_rows(doc, tissue, spec):
    fam = "primary_family" if spec == "full" else "sensitivity_family"
    return doc["pooled"][tissue][fam]["rows"]


def frozen_analysis(doc):
    out = {}
    for tissue in ("brain", "pbmc"):
        for spec in ("full", "non_degree"):
            rows = family_rows(doc, tissue, spec)
            labels = [r["model_label"] for r in rows]
            p = {nt: np.array([r[key]["p_mc"] for r in rows]) for nt, key in NULL_KEYS.items()}
            stored_q = {nt: np.array([r[key]["bh_q_family"] for r in rows]) for nt, key in NULL_KEYS.items()}
            res = {"labels": labels,
                   "rho": [r["observed_partial_rho"] for r in rows],
                   "p": {nt: p[nt].tolist() for nt in p}, "q": {}, "support": {}}
            for method in ("BH", "BY", "Holm"):
                qm, qd = adjust(p["mantel"], method), adjust(p["degree"], method)
                if method == "BH":
                    for nt, q in (("mantel", qm), ("degree", qd)):
                        assert np.allclose(np.round(q, 6), stored_q[nt], atol=1e-6), (tissue, spec, nt)
                res["q"][method] = {"mantel": qm.tolist(), "degree": qd.tolist()}
                res["support"][method] = ((qm < ALPHA) & (qd < ALPHA)).tolist()
            iut = np.maximum(p["mantel"], p["degree"])
            for method in ("BH", "BY"):
                q = adjust(iut, method)
                res["q"][f"IUT-{method}"] = q.tolist()
                res["support"][f"IUT-{method}"] = (q < ALPHA).tolist()
            out[f"{tissue}_{spec}"] = res
    return out


def joint_null(res_dir, tissue, spec, null_type, labels):
    """(N, m) matrix of null statistics in family row order."""
    cols = []
    groups = ["brain"] if tissue == "brain" else ["pbmc_shared", "pbmc_scgpt", "pbmc_uce"]
    arrays = {}
    for g in groups:
        f = res_dir / f"fm_{g}_{null_type}_{spec}.npz"
        Z = np.load(f)
        for k in Z.files:
            if not k.startswith("_"):
                arrays[k] = Z[k]
    for lab in labels:
        cols.append(arrays[lab])
    return np.column_stack(cols)


def deep_analysis(doc, res_dir):
    obs_exact = common.observed_unrounded(common.load_groups())
    out = {}
    for tissue in ("brain", "pbmc"):
        for spec in ("full", "non_degree"):
            rows = family_rows(doc, tissue, spec)
            labels = [r["model_label"] for r in rows]
            obs = np.array([obs_exact[(tissue, spec, lab)] for lab in labels])
            res = {"labels": labels, "p": {}, "q": {}, "support": {}}
            p = {}
            for nt in NULL_KEYS:
                M = joint_null(res_dir, tissue, spec, nt, labels)
                raw, wy = wy_stepdown_minp(M, obs)
                p[nt] = raw
                res["p"][nt] = raw.tolist()
                res["q"].setdefault("WY", {})[nt] = wy.tolist()
                res["n_perm"] = int(M.shape[0])
                res.setdefault("null_sd", {})[nt] = M.std(axis=0).tolist()
            res["support"]["WY"] = ((np.array(res["q"]["WY"]["mantel"]) < ALPHA)
                                    & (np.array(res["q"]["WY"]["degree"]) < ALPHA)).tolist()
            for method in ("BH", "BY", "Holm"):
                qm, qd = adjust(p["mantel"], method), adjust(p["degree"], method)
                res["q"][method] = {"mantel": qm.tolist(), "degree": qd.tolist()}
                res["support"][method] = ((qm < ALPHA) & (qd < ALPHA)).tolist()
            iut = np.maximum(p["mantel"], p["degree"])
            for method in ("BH", "BY"):
                q = adjust(iut, method)
                res["q"][f"IUT-{method}"] = q.tolist()
                res["support"][f"IUT-{method}"] = (q < ALPHA).tolist()
            out[f"{tissue}_{spec}"] = res
    return out


def counts(block):
    methods = sorted({m for fam in block.values() for m in fam["support"]})
    return {m: {spec: int(sum(sum(block[f"{t}_{spec}"]["support"][m]) for t in ("brain", "pbmc")))
                for spec in ("full", "non_degree")} for m in methods}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary-only", action="store_true",
                    help="Check frozen p-values only; do not load graph or null-array caches")
    ap.add_argument("--output", type=Path,
                    default=Path(common.fpa.ROOT) / "output/revision_checks/multiplicity_robustness.json")
    args = ap.parse_args()
    doc = json.loads(common.AUDIT_JSON.read_text())
    frozen = frozen_analysis(doc)
    result = {
        "schema_version": 1,
        "analysis": "multiplicity_robustness",
        "alpha": ALPHA,
        "rule_note": ("Dual-null Support under a method = both nulls below alpha after that "
                      "method within the frozen family; IUT = one test per row on max(pM, pD)."),
        "frozen_n999": {"families": frozen, "support_counts": counts(frozen)},
        "inputs": {"fixed_panel_audit_v2_sha256": common.sha256_file(common.AUDIT_JSON)},
    }
    res_dir = common.OUT_DIR / "resolution_n9999"
    if not args.summary_only and all((res_dir / f"fm_{g}_{nt}_{sp}.npz").exists()
           for g in ("brain", "pbmc_shared", "pbmc_scgpt", "pbmc_uce")
           for nt in ("mantel", "degree") for sp in ("full", "non_degree")):
        deep = deep_analysis(doc, res_dir)
        result["resolution_n9999"] = {"families": deep, "support_counts": counts(deep)}
    else:
        log("Frozen summary only; N=9,999 graph-based analysis was not executed")
    common.write_json(args.output, result)
    log("support counts (frozen N=999):", result["frozen_n999"]["support_counts"])
    if "resolution_n9999" in result:
        log("support counts (N=9,999):", result["resolution_n9999"]["support_counts"])


if __name__ == "__main__":
    main()
