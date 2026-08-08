#!/usr/bin/env python
"""Task 2 of the reg-audit CPU must-fixes: the specification gate.

Two questions, two JSONs, no GPU.

(1) marginal_vs_adjusted_v2.json -- the nested covariate ladder for every
    (tissue, model) on the fixed panel's OWN edge set and with the fixed panel's
    OWN statistic, so the rungs are comparable to the headline numbers rather
    than merely similar to them. confound_ablation_v2.json answered a narrower
    version of this on the pre-marker-mask edge set for 5 of the 9 model rows;
    it is reloaded here and cross-checked, not replaced.

(2) spec_sensitivity_v2.json -- the primary statistic recomputed with and
    without the covariates that are derived from the outcome graph itself
    (tf_outdeg, atac_indeg, both functions of G_ATAC). Conditioning an estimate
    on functions of its own outcome is the failure mode that gates the GPU
    spend: if the headline effects are an artifact of that conditioning,
    replicating the statistic on more model families multiplies the artifact.

Rungs (an FM row, use_coexp=True):
    marginal              no controls
    coexp_only            + rank(co-expression)
    nondegree_only        + peakcount, genelen, detection, GC   (no coexp)
    degree_only           + tf_outdeg, atac_indeg               (no coexp)
    coexp_plus_nondegree  + coexp + non-degree      == audit spec 'non_degree'
    coexp_plus_full       + coexp + non-degree + degree  == audit spec 'full'
The two isolating rungs (nondegree_only, degree_only) exist to attribute the
marginal->full drop to a covariate block rather than merely observing it.

For the co-expression baseline row the coexp rungs are skipped: partialling
co-expression out of co-expression is degenerate by construction.

Monte-Carlo p-values are NOT recomputed here -- they are read from
fixed_panel_audit_v2.json, whose permutation nulls were drawn under these exact
two specs. CPU-only, reuses cached graphs.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixed_panel_audit as fpa  # noqa: E402
import run_fixed_panel_audit as rfp  # noqa: E402

OUT = fpa.OUT
UTC = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

# Covariates computed from G_ATAC, i.e. from the outcome of the comparison.
OUTCOME_DERIVED = ["tf_outdeg", "atac_indeg"]


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def _controls(spec: str, ii, jj, peakcount, genelen, detv, gc, tf_outdeg, atac_indeg,
              co_rank: Optional[np.ndarray]) -> Optional[np.ndarray]:
    """Control block for one rung. None means 'no controls' (plain correlation)."""
    cols: List[np.ndarray] = []
    if co_rank is not None:
        cols.append(co_rank)
    if spec in ("nondegree", "full"):
        cols += [fpa.zscore(peakcount[jj]), fpa.zscore(genelen[jj]),
                 fpa.zscore(detv[jj]), fpa.zscore(gc[jj])]
    if spec in ("degree", "full"):
        cols += [fpa.zscore(tf_outdeg[ii]), fpa.zscore(atac_indeg[jj])]
    if not cols:
        return None
    return np.column_stack(cols)


def ladder(fm_v, atac_v, co_v, ii, jj, peakcount, genelen, detv, gc,
           tf_outdeg, atac_indeg, use_coexp: bool) -> Dict[str, Optional[float]]:
    """The nested ladder for one (tissue, model) row."""
    r_fm, r_atac, r_co = rankdata(fm_v), rankdata(atac_v), rankdata(co_v)

    def rung(spec: str, with_coexp: bool) -> Optional[float]:
        if with_coexp and not use_coexp:
            return None  # degenerate: fm IS co-expression here
        C = _controls(spec, ii, jj, peakcount, genelen, detv, gc,
                      tf_outdeg, atac_indeg, r_co if with_coexp else None)
        if C is None:
            return float(np.corrcoef(r_fm, r_atac)[0, 1])
        return fpa.pcorr(r_fm, r_atac, C)

    return {
        "marginal": rung("none", False),
        "coexp_only": rung("none", True),
        "nondegree_only": rung("nondegree", False),
        "degree_only": rung("degree", False),
        "coexp_plus_nondegree": rung("nondegree", True),
        "coexp_plus_full": rung("full", True),
    }


def _round(v: Optional[float], nd: int = 6) -> Optional[float]:
    return None if v is None else round(float(v), nd)


def _sign(v: Optional[float]) -> Optional[int]:
    return None if v is None else int(np.sign(v))


def run_tissue(tissue: str, atac_file: str, G_atac, co, models: Dict[str, np.ndarray],
               tf_rows: np.ndarray) -> List[dict]:
    Ng = G_atac.shape[0]
    peakcount, genelen, gc, detv = fpa.build_confounds(atac_file)
    tf_outdeg = (G_atac > 0).sum(1).astype(np.float32)
    atac_indeg = (G_atac > 0).sum(0).astype(np.float32)

    genes = json.loads(Path(fpa.MANI).read_text())["genes"]
    ii_all = np.repeat(tf_rows, Ng)
    jj_all = np.tile(np.arange(Ng), len(tf_rows))
    m = fpa.edge_mask(tissue, genes, tf_rows, ii_all, jj_all)
    ii, jj = ii_all[m], jj_all[m]

    atac_v = G_atac[ii, jj]
    co_v = co[ii, jj]
    log(f"  {tissue}: {len(ii)} edges, {len(tf_rows)} TF rows, marker_mask="
        f"{tissue in fpa.MARKER_TISSUES}")

    rows = []
    entries: List[tuple] = [("co_expression", co, False)]
    entries += [(label, G, True) for label, G in models.items()]

    for label, G_fm, use_coexp in entries:
        fm_v = G_fm[ii, jj]
        st = ladder(fm_v, atac_v, co_v, ii, jj, peakcount, genelen, detv, gc,
                    tf_outdeg, atac_indeg, use_coexp)
        primary = st["coexp_plus_full"] if use_coexp else st["degree_only"]
        marg = st["marginal"]
        # Both rungs are unconditionally computed for every row shape reaching here:
        # 'marginal' takes no controls, and 'primary' is chosen to match use_coexp.
        assert primary is not None and marg is not None, f"empty primary rung for {label}"
        rows.append({
            "tissue": tissue,
            "model_label": label,
            "is_baseline": not use_coexp,
            "n_pairs": int(len(ii)),
            **{k: _round(v) for k, v in st.items()},
            "delta_marginal_to_primary": _round(primary - marg),
            "sign_flip_marginal_to_primary": bool(_sign(marg) != _sign(primary)),
            "primary_rung": "coexp_plus_full" if use_coexp else "degree_only",
        })
        log(f"    [{label}] marg={marg:+.4f}"
            f" coexp={st['coexp_only'] if st['coexp_only'] is None else format(st['coexp_only'], '+.4f')}"
            f" nd_only={st['nondegree_only']:+.4f} deg_only={st['degree_only']:+.4f}"
            f" -> nd={st['coexp_plus_nondegree'] if st['coexp_plus_nondegree'] is None else format(st['coexp_plus_nondegree'], '+.4f')}"
            f" full={st['coexp_plus_full'] if st['coexp_plus_full'] is None else format(st['coexp_plus_full'], '+.4f')}")
    return rows


def cross_check_ablation(rows: List[dict]) -> dict:
    """confound_ablation_v2.json ran the same idea on a DIFFERENT edge set (no
    marker mask) for 5 brain rows. Agreement is expected to be close but not
    exact; the point of the check is that nothing moved by more than rounding
    plus the 446-edge mask, and that no rung changed sign."""
    path = f"{OUT}/confound_ablation_v2.json"
    if not os.path.exists(path):
        return {"status": "absent", "path": path}
    old = json.loads(Path(path).read_text())
    alias = {"coexp_vs_atac": "co_expression", "geneformer_embed": "geneformer_embed",
             "scfoundation": "scFoundation_encoder", "uce": "UCE_encoder",
             "geneformer_ko_raw": "geneformer_ko_raw"}
    by = {(r["tissue"], r["model_label"]): r for r in rows}
    comps = []
    for o in old:
        label = alias.get(o["label"])
        new = by.get(("brain", label))
        if new is None:
            comps.append({"ablation_label": o["label"], "status": "no_matching_row"})
            continue
        deltas = {}
        for rung in ("marginal", "coexp_only", "coexp_plus_nondegree", "coexp_plus_full"):
            ov, nv = o.get(rung), new.get(rung)
            if ov is None or nv is None:
                deltas[rung] = None
                continue
            deltas[rung] = {
                "ablation": ov, "recomputed": nv,
                "abs_delta": round(abs(nv - ov), 6),
                "sign_agrees": bool(np.sign(ov) == np.sign(nv)),
            }
        finite = [d for d in deltas.values() if d is not None]
        comps.append({
            "ablation_label": o["label"], "model_label": label, "status": "compared",
            "rungs": deltas,
            "max_abs_delta": round(max(d["abs_delta"] for d in finite), 6) if finite else None,
            "all_signs_agree": bool(all(d["sign_agrees"] for d in finite)),
        })
    done = [c for c in comps if c["status"] == "compared"]
    return {
        "status": "compared",
        "note": ("confound_ablation_v2.json used the unmasked brain edge set (534754 pairs); "
                 "this recomputation uses the fixed panel's marker-masked set (534308). "
                 "Deltas are expected to be small but non-zero."),
        "n_compared": len(done),
        "max_abs_delta_across_rows": round(max((c["max_abs_delta"] for c in done
                                                if c["max_abs_delta"] is not None), default=0.0), 6),
        "all_signs_agree": bool(all(c["all_signs_agree"] for c in done)),
        "comparisons": comps,
    }


def _audit_pvalues() -> Dict[tuple, dict]:
    """Mantel + degree-preserving p-values from the fixed-panel audit, keyed by
    (tissue, model_label, confound_spec). Not recomputed here."""
    path = f"{OUT}/fixed_panel_audit_v2.json"
    if not os.path.exists(path):
        return {}
    doc = json.loads(Path(path).read_text())
    out = {}
    for tissue, blk in doc.get("pooled", {}).items():
        for r in blk.get("rows", []):
            if r["row_type"] != "pooled_fm":
                continue
            out[(tissue, r["model_label"], r["confound_spec"])] = {
                "observed_partial_rho": r.get("observed_partial_rho"),
                "mantel_p": (r.get("mantel") or {}).get("p_mc"),
                "degree_preserving_p": (r.get("degree_preserving") or {}).get("p_mc"),
            }
    return out


def spec_sensitivity(rows: List[dict]) -> dict:
    """full (conditions on outcome-derived degree) vs non_degree (does not)."""
    pv = _audit_pvalues()
    fm = [r for r in rows if not r["is_baseline"]]
    base = {r["tissue"]: r for r in rows if r["is_baseline"]}

    out_rows = []
    for r in fm:
        t, label = r["tissue"], r["model_label"]
        full, nd = r["coexp_plus_full"], r["coexp_plus_nondegree"]
        a_full = pv.get((t, label, "full"), {})
        a_nd = pv.get((t, label, "non_degree"), {})
        coexp_marg = base[t]["marginal"]
        out_rows.append({
            "tissue": t,
            "model_label": label,
            "rho_full": full,
            "rho_non_degree": nd,
            "delta_full_minus_non_degree": _round(full - nd),
            "sign_flip_between_specs": bool(_sign(full) != _sign(nd)),
            "mantel_p_full": a_full.get("mantel_p"),
            "mantel_p_non_degree": a_nd.get("mantel_p"),
            "degree_preserving_p_full": a_full.get("degree_preserving_p"),
            "degree_preserving_p_non_degree": a_nd.get("degree_preserving_p"),
            "audit_rho_full_matches": (
                None if a_full.get("observed_partial_rho") is None
                else bool(abs(a_full["observed_partial_rho"] - full) < 1e-6)),
            "audit_rho_non_degree_matches": (
                None if a_nd.get("observed_partial_rho") is None
                else bool(abs(a_nd["observed_partial_rho"] - nd) < 1e-6)),
            "abs_rho_full_vs_coexp_marginal": _round(abs(full) / abs(coexp_marg)),
            "abs_rho_non_degree_vs_coexp_marginal": _round(abs(nd) / abs(coexp_marg)),
        })

    # Does the spec choice reorder the models within a tissue?
    rank_stability = {}
    for t in sorted({r["tissue"] for r in out_rows}):
        sub = [r for r in out_rows if r["tissue"] == t]
        if len(sub) < 3:
            rank_stability[t] = {"n_models": len(sub),
                                 "spearman_full_vs_non_degree": None,
                                 "note": "fewer than 3 models; rank correlation not reported"}
            continue
        f = [r["rho_full"] for r in sub]
        n = [r["rho_non_degree"] for r in sub]
        rank_stability[t] = {
            "n_models": len(sub),
            "spearman_full_vs_non_degree": round(float(spearmanr(f, n).statistic), 4),
            "argmax_full": max(sub, key=lambda r: r["rho_full"])["model_label"],
            "argmax_non_degree": max(sub, key=lambda r: r["rho_non_degree"])["model_label"],
        }

    # Significance count under each spec, and whether the two null models agree.
    # They are different nulls: Mantel permutes gene labels, the degree-preserving
    # null shuffles within degree strata. Disagreement is informative, not a bug.
    def _sig(key: str) -> List[str]:
        return [f"{r['tissue']}/{r['model_label']}" for r in out_rows
                if r.get(key) is not None and r[key] < 0.05]

    sig = {
        "mantel_full": _sig("mantel_p_full"),
        "mantel_non_degree": _sig("mantel_p_non_degree"),
        "degree_preserving_full": _sig("degree_preserving_p_full"),
        "degree_preserving_non_degree": _sig("degree_preserving_p_non_degree"),
    }
    null_agreement = {
        "note": ("Mantel permutes gene labels; the degree-preserving null shuffles within "
                 "degree strata. Under 'full' they agree. Under 'non_degree' they do not: "
                 "Mantel finds nothing while the degree-preserving null still flags rows. "
                 "With no degree covariate in the design, the degree-preserving null no "
                 "longer breaks the degree structure the statistic is picking up, so the "
                 "more conservative Mantel result governs."),
        "n_significant_at_0p05": {k: len(v) for k, v in sig.items()},
        "significant_rows": sig,
        "nulls_agree_under_full": sorted(sig["mantel_full"]) == sorted(sig["degree_preserving_full"]),
        "nulls_agree_under_non_degree": (
            sorted(sig["mantel_non_degree"]) == sorted(sig["degree_preserving_non_degree"])),
    }

    flips = [r for r in out_rows if r["sign_flip_between_specs"]]
    exceeds = [r for r in out_rows
               if max(abs(r["rho_full"]), abs(r["rho_non_degree"]))
               >= abs(base[r["tissue"]]["marginal"])]
    return {
        "question": ("Does the primary statistic depend on conditioning on covariates "
                     "derived from the outcome graph G_ATAC?"),
        "outcome_derived_covariates": OUTCOME_DERIVED,
        "specs": {
            "full": "coexp + peakcount, genelen, detection, GC + tf_outdeg, atac_indeg (outcome-derived)",
            "non_degree": "coexp + peakcount, genelen, detection, GC (no outcome-derived covariate)",
        },
        "rows": out_rows,
        "rank_stability": rank_stability,
        "null_model_agreement": null_agreement,
        "summary": {
            "n_model_rows": len(out_rows),
            "n_sign_flips_between_specs": len(flips),
            "sign_flipping_rows": [f"{r['tissue']}/{r['model_label']}" for r in flips],
            "max_abs_delta_between_specs": round(
                max(abs(r["delta_full_minus_non_degree"]) for r in out_rows), 6),
            "n_rows_reaching_coexp_marginal_under_either_spec": len(exceeds),
            "coexp_marginal_by_tissue": {t: base[t]["marginal"] for t in sorted(base)},
            "all_audit_rho_reproduced": bool(all(
                (r["audit_rho_full_matches"] in (True, None))
                and (r["audit_rho_non_degree_matches"] in (True, None)) for r in out_rows)),
        },
    }


def provenance(paths: Dict[str, str]) -> dict:
    return {
        "generated_utc": UTC,
        "manifest_sha256": fpa.load_manifest()[2],
        "inputs": {k: {"path": v, "sha256": fpa.sha256_file(v)}
                   for k, v in paths.items() if os.path.exists(v)},
        "code": {"module": "src/v2/spec_audit_v2.py",
                 "reused": ["fixed_panel_audit.partial_rho ladder via fpa.pcorr",
                            "fixed_panel_audit.edge_mask",
                            "fixed_panel_audit.build_confounds"]},
    }


def main() -> None:
    ATAC_B = f"{rfp.DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
    ATAC_P = f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"
    for p in (ATAC_B, ATAC_P, fpa.MANI, fpa.COORDS, fpa.HG38):
        if not os.path.exists(p):
            raise FileNotFoundError(f"required input missing: {p}")

    log("=== loading cached graphs ===")
    G_atac_b, co_b, models_b, tf_b, _ = rfp.load_pooled_brain()
    G_atac_p, co_p, models_p, tf_p, _ = rfp.load_pooled_pbmc()
    ko_models, _ = rfp.load_geneformer_ko()
    models_b = {**models_b, **ko_models}

    log("=== ladder: marginal -> adjusted ===")
    rows = run_tissue("brain", ATAC_B, G_atac_b, co_b, models_b, tf_b)
    rows += run_tissue("pbmc", ATAC_P, G_atac_p, co_p, models_p, tf_p)

    prov = provenance({
        "atac_brain": ATAC_B, "atac_pbmc": ATAC_P,
        "fmgraphs_pooled_brain": f"{OUT}/fmgraphs_pooled_v2.npz",
        "fmgraphs_pooled_pbmc": f"{OUT}/pbmc_fmgraphs_pooled.npz",
        "G_atac_brain": f"{OUT}/G_ATAC_v2_GSE174367.npz",
        "G_atac_pbmc": f"{OUT}/G_ATAC_v2_PBMC10k.npz",
        "G_ko": f"{OUT}/G_ko_v2.npz",
    })

    doc1 = {
        "schema_version": 1,
        "generated_utc": UTC,
        "purpose": ("Nested covariate ladder for every (tissue, model) on the fixed panel's "
                    "edge set, using the fixed panel's statistic."),
        "rungs": {
            "marginal": "no controls",
            "coexp_only": "+ rank(co-expression)",
            "nondegree_only": "+ peakcount, genelen, detection, GC (no coexp)",
            "degree_only": "+ tf_outdeg, atac_indeg (no coexp) -- both outcome-derived",
            "coexp_plus_nondegree": "audit spec 'non_degree'",
            "coexp_plus_full": "audit spec 'full' (headline)",
        },
        "baseline_note": ("The co_expression row is the baseline, not a model; its coexp rungs "
                          "are null because partialling co-expression out of itself is degenerate. "
                          "Its primary rung is degree_only."),
        "rows": rows,
        "cross_check_confound_ablation": cross_check_ablation(rows),
        "provenance": prov,
    }
    Path(f"{OUT}/marginal_vs_adjusted_v2.json").write_text(json.dumps(doc1, indent=2))
    log(f"SAVED {OUT}/marginal_vs_adjusted_v2.json")

    doc2 = {"schema_version": 1, "generated_utc": UTC, **spec_sensitivity(rows),
            "provenance": prov}
    Path(f"{OUT}/spec_sensitivity_v2.json").write_text(json.dumps(doc2, indent=2))
    log(f"SAVED {OUT}/spec_sensitivity_v2.json")

    s = doc2["summary"]
    log(f"sign flips between specs: {s['n_sign_flips_between_specs']}/{s['n_model_rows']}"
        f" {s['sign_flipping_rows']}")
    log(f"max |delta| between specs: {s['max_abs_delta_between_specs']}")
    log(f"rows reaching the co-expression marginal under either spec: "
        f"{s['n_rows_reaching_coexp_marginal_under_either_spec']}")
    log(f"audit rho reproduced: {s['all_audit_rho_reproduced']}")


if __name__ == "__main__":
    main()
