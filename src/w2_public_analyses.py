#!/usr/bin/env python3
"""Wave-2 analyses that run on *public* audit JSON only (capsule-safe).

Produces:
  results/fm_vs_baseline_observed_v2.public.json   (0x21)
  results/dual_null_oc_independence_v2.public.json (0x22)
  results/tf_probe_contrasts_no_floor_v2.public.json (0x25)
  results/nondegree_null_pattern_v2.public.json    (0x23 support table)

No NPZ / H5AD / model weights required.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RNG = np.random.default_rng(20260809)


def bh_q(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=pvalues.__getitem__)
    out = [1.0] * n
    running = 1.0
    for reverse_index in range(n - 1, -1, -1):
        index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, pvalues[index] * n / rank, 1.0)
        out[index] = running
    return out


def load_json(name: str):
    return json.loads((RES / name).read_text())


def fm_vs_baseline() -> dict:
    audit = load_json("fixed_panel_audit_v2.public.json")
    bases = {
        "brain": load_json("brain_coexp_baseline_null_v2.public.json"),
        "pbmc": load_json("pbmc_coexp_baseline_null_v2.public.json"),
    }
    rows = []
    for tissue, base in bases.items():
        br = float(base["observed_rho"])
        for r in audit["pooled"][tissue]["rows"]:
            if r.get("row_type") != "pooled_fm" or r.get("confound_spec") != "full":
                continue
            rho = float(r["observed_partial_rho"])
            dual = (
                float(r["mantel"]["bh_q_family"]) < 0.05
                and float(r["degree_preserving"]["bh_q_family"]) < 0.05
            )
            delta = rho - br
            rows.append(
                {
                    "tissue": tissue,
                    "model_label": r["model_label"],
                    "rho_full": rho,
                    "rho_baseline": br,
                    "delta_rho": delta,
                    "beats_baseline": bool(delta > 0),
                    "dual_null_full": dual,
                    "dual_and_beats_baseline": bool(dual and delta > 0),
                }
            )
    dual_n = sum(1 for x in rows if x["dual_null_full"])
    dual_beats = sum(1 for x in rows if x["dual_and_beats_baseline"])
    beats = sum(1 for x in rows if x["beats_baseline"])
    return {
        "schema_version": 1,
        "method": "observed_delta_v1",
        "note": (
            "Observed FM−baseline partial-ρ differences on the fixed panel. "
            "Shared-null randomization of Δ requires monorepo NPZ caches; not run here."
        ),
        "n_fm_rows": len(rows),
        "n_beats_baseline": beats,
        "n_dual_null_full": dual_n,
        "n_dual_and_beats_baseline": dual_beats,
        "summary": (
            f"{dual_beats}/{dual_n} dual-null full rows have ρ_full > tissue co-expression baseline; "
            f"{beats}/{len(rows)} FM rows beat baseline ignoring dual-null."
        ),
        "rows": rows,
    }


def dual_null_oc(n_sims: int = 20000) -> dict:
    """Independence reference OC for dual-null Support under global null."""
    audit = load_json("fixed_panel_audit_v2.public.json")
    family_sizes = {}
    for tissue in ("brain", "pbmc"):
        n = sum(
            1
            for r in audit["pooled"][tissue]["rows"]
            if r.get("row_type") == "pooled_fm" and r.get("confound_spec") == "full"
        )
        family_sizes[tissue] = n

    out_families = {}
    for tissue, n in family_sizes.items():
        # Per simulation: n independent p_M, n independent p_D; BH within null family.
        any_dual = 0
        mean_dual_rows = 0.0
        for _ in range(n_sims):
            p_m = RNG.random(n)
            p_d = RNG.random(n)
            q_m = np.array(bh_q(p_m.tolist()))
            q_d = np.array(bh_q(p_d.tolist()))
            dual = (q_m < 0.05) & (q_d < 0.05)
            k = int(dual.sum())
            mean_dual_rows += k
            if k > 0:
                any_dual += 1
        out_families[tissue] = {
            "n_rows_in_family": n,
            "n_sims": n_sims,
            "rate_at_least_one_dual_null_row": any_dual / n_sims,
            "mean_dual_null_rows_per_sim": mean_dual_rows / n_sims,
            "assumption": "independent Uniform p_M, p_D; separate BH per null; no shared-batch dependence",
        }

    # Observed dual-null counts under full
    observed = {}
    for tissue in ("brain", "pbmc"):
        k = 0
        for r in audit["pooled"][tissue]["rows"]:
            if r.get("row_type") != "pooled_fm" or r.get("confound_spec") != "full":
                continue
            if (
                float(r["mantel"]["bh_q_family"]) < 0.05
                and float(r["degree_preserving"]["bh_q_family"]) < 0.05
            ):
                k += 1
        observed[tissue] = k

    return {
        "schema_version": 1,
        "method": "global_null_independence_mc_v1",
        "note": (
            "Reference operating characteristic under independent null p-values. "
            "The audit shares one proxy null batch across rows, inducing positive dependence; "
            "this OC is not exact FDR under the audit stream."
        ),
        "observed_dual_null_full_counts": observed,
        "families": out_families,
        "interpretation": (
            "If independence OC rate of ≥1 dual-null row is far below the event of observing "
            f"{sum(observed.values())} dual-null rows across tissues, dual-null Support is not "
            "explained by naive independent double-dipping alone — still not causal recovery."
        ),
    }


def nondegree_null_pattern() -> dict:
    audit = load_json("fixed_panel_audit_v2.public.json")
    rows = []
    counts = {"dual": 0, "M_only": 0, "D_only": 0, "neither": 0}
    for tissue in ("brain", "pbmc"):
        for r in audit["pooled"][tissue]["rows"]:
            if r.get("row_type") != "pooled_fm" or r.get("confound_spec") != "non_degree":
                continue
            m = float(r["mantel"]["bh_q_family"]) < 0.05
            d = float(r["degree_preserving"]["bh_q_family"]) < 0.05
            if m and d:
                tag = "dual"
            elif m and not d:
                tag = "M_only"
            elif d and not m:
                tag = "D_only"
            else:
                tag = "neither"
            counts[tag] += 1
            rows.append(
                {
                    "tissue": tissue,
                    "model_label": r["model_label"],
                    "rho": float(r["observed_partial_rho"]),
                    "q_M": float(r["mantel"]["bh_q_family"]),
                    "q_D": float(r["degree_preserving"]["bh_q_family"]),
                    "pattern": tag,
                }
            )
    return {
        "schema_version": 1,
        "confound_spec": "non_degree",
        "counts": counts,
        "rows": rows,
        "note": (
            "Under non-degree, all rejections in audit v2 are D-only (row-shuffle). "
            "Dual-Support never fires; D-only must not be reported as dual-null Support."
        ),
    }


def probe_no_floor() -> dict:
    probe = load_json("tf_probe_pair_stats_v2.public.json")
    contrasts = probe["contrasts_vs_baseline"]
    # Exclude random_floor from the BH family for FM vs coexp contrasts.
    keys = [k for k in contrasts if k != "random_floor"]
    pvals = [float(contrasts[k]["signflip_p"]) for k in keys]
    qs = bh_q(pvals)
    out_rows = {}
    for k, q, p in zip(keys, qs, pvals):
        c = dict(contrasts[k])
        c["signflip_q_family_without_random_floor"] = q
        c["significant_q05_without_random_floor"] = bool(q < 0.05)
        c["signflip_p"] = p
        out_rows[k] = c
    # Also report original random_floor for transparency
    floor = dict(contrasts["random_floor"])
    floor["note"] = "Excluded from BH family in sensitivity; original q retained for reference"
    return {
        "schema_version": 1,
        "method": "bh_within_contrast_family_excluding_random_floor_v1",
        "n_perm": probe.get("n_perm"),
        "tissue": "PBMC only (brain not probed)",
        "contrasts_without_random_floor": out_rows,
        "random_floor_original": floor,
        "summary": {
            k: {
                "signflip_q_no_floor": out_rows[k]["signflip_q_family_without_random_floor"],
                "sig": out_rows[k]["significant_q05_without_random_floor"],
                "original_signflip_q": float(contrasts[k]["signflip_q"]),
            }
            for k in keys
        },
    }


def main():
    RES.mkdir(parents=True, exist_ok=True)
    artifacts = {
        "fm_vs_baseline_observed_v2.public.json": fm_vs_baseline(),
        "dual_null_oc_independence_v2.public.json": dual_null_oc(),
        "nondegree_null_pattern_v2.public.json": nondegree_null_pattern(),
        "tf_probe_contrasts_no_floor_v2.public.json": probe_no_floor(),
    }
    for name, obj in artifacts.items():
        path = RES / name
        path.write_text(json.dumps(obj, indent=2) + "\n")
        print("wrote", path.relative_to(ROOT), "bytes", path.stat().st_size)
    # print short summaries
    print("0x21", artifacts["fm_vs_baseline_observed_v2.public.json"]["summary"])
    for t, fam in artifacts["dual_null_oc_independence_v2.public.json"]["families"].items():
        print(
            "0x22",
            t,
            "P(≥1 dual)≈",
            round(fam["rate_at_least_one_dual_null_row"], 4),
            "obs dual",
            artifacts["dual_null_oc_independence_v2.public.json"]["observed_dual_null_full_counts"][t],
        )
    print("0x23 counts", artifacts["nondegree_null_pattern_v2.public.json"]["counts"])
    print("0x25", artifacts["tf_probe_contrasts_no_floor_v2.public.json"]["summary"])


if __name__ == "__main__":
    main()
