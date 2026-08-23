#!/usr/bin/env python3
"""Construct-lane hooks: Mantel / additive-decomp against locked proxies.

Runs on **local** ``G_ATAC`` NPZ when present (fibroblast / brain / PBMC freeze
caches). Does **not** download data, does **not** mutate PeerJ Support public
JSON or MANIFEST locks. Heavy outputs land under ``results/v2/extension/``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from paths import (  # noqa: E402
    EXTENSION_ROOT,
    HEAVY_ARTIFACT_ROOT,
    PLAN_ROOT,
    ROOT,
    local_asset_report,
    redact_path,
    resolve_g_atac_npz,
    resolve_local_atac,
)
from registry import load_extension_registry  # noqa: E402

COMPARE_TAGS = ("GSE174367", "PBMC10k", "GSE206767")
N_ADDITIVE_ITER = 50


def fibro_style_env(tissue_id: str, atac_file: str | None = None) -> dict[str, str]:
    """Return env vars mirroring the GSE206767 construct seam."""
    reg = load_extension_registry()
    reg.assert_may_emit_g_atac(tissue_id)
    meta = reg.get_tissue(tissue_id)
    tag = str(meta.get("g_atac_tag") or tissue_id)
    resolved = resolve_local_atac(tissue_id, explicit=atac_file)
    atac_redacted = redact_path(resolved) if resolved else redact_path(atac_file) or ""
    env = {
        "TAG": tag,
        "ATAC_FILE": atac_redacted,
        "SCREG_EXTENSION_OUT": f"results/v2/extension/construct/{tag}",
        "SCREG_EXTENSION_LANE": "construct",
        "SCREG_PEERJ_SUPPORT_LOCK": "1",
    }
    return env


def _load_consensus(path: Path) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(path, allow_pickle=False)
    types = [str(t) for t in z["types"]]
    g = np.mean([z[f"G_{t}"] for t in types], axis=0).astype(np.float64)
    return g, np.asarray(z["tf_rows"], dtype=int)


def _additive_fit(m: np.ndarray, tf_rows: np.ndarray, n_iter: int = N_ADDITIVE_ITER):
    sub = m[tf_rows, :].copy()
    mask = np.ones_like(sub, dtype=bool)
    for k, tf in enumerate(tf_rows):
        mask[k, int(tf)] = False
    mu = float(sub[mask].mean())
    r = np.zeros(len(tf_rows))
    c = np.zeros(m.shape[1])
    for _ in range(n_iter):
        r_new = np.where(
            mask.sum(1) > 0,
            (sub - mu - c[None, :]).sum(1) / mask.sum(1),
            0.0,
        )
        c_new = np.where(
            mask.sum(0) > 0,
            (sub - mu - r_new[:, None]).sum(0) / mask.sum(0),
            0.0,
        )
        if np.allclose(r_new, r) and np.allclose(c_new, c):
            r, c = r_new, c_new
            break
        r, c = r_new, c_new
    return mu, r, c


def _pair_mantel_decomp(
    tag_a: str,
    ga: np.ndarray,
    tf_a: np.ndarray,
    tag_b: str,
    gb: np.ndarray,
    tf_b: np.ndarray,
) -> dict[str, Any]:
    tf_common = np.intersect1d(tf_a, tf_b)
    ng = ga.shape[0]
    ii = np.repeat(tf_common, ng)
    jj = np.tile(np.arange(ng), len(tf_common))
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    x, y = ga[ii, jj], gb[ii, jj]
    observed = float(spearmanr(x, y).statistic)
    mu, r, c = _additive_fit(ga, tf_common)
    row_of = {int(tf): k for k, tf in enumerate(tf_common)}
    pred = mu + r[[row_of[int(t)] for t in ii]] + c[jj]
    additive_pred_rho = float(spearmanr(pred, y).statistic)
    mu_b, r_b, c_b = _additive_fit(gb, tf_common)
    pred_b = mu_b + r_b[[row_of[int(t)] for t in ii]] + c_b[jj]
    resid_rho = float(spearmanr(x - pred, y - pred_b).statistic)
    return {
        "pair": [tag_a, tag_b],
        "n_tf_common": int(len(tf_common)),
        "observed_spearman": round(observed, 6),
        "additive_pred_spearman": round(additive_pred_rho, 6),
        "residual_spearman_after_own_additive_fits": round(resid_rho, 6),
        "fraction_explained_by_additive_marginals": (
            round(additive_pred_rho / observed, 4) if observed else None
        ),
    }


def run_construct(
    tissue_id: str,
    *,
    atac_file: str | None = None,
    execute: bool = False,
) -> dict[str, Any]:
    """Plan and optionally execute Mantel/decomp for a construct-lane tissue."""
    reg = load_extension_registry()
    reg.assert_may_emit_g_atac(tissue_id)
    meta = reg.get_tissue(tissue_id)
    tag = str(meta.get("g_atac_tag") or tissue_id)
    assets = local_asset_report(tissue_id, tag)
    if atac_file:
        assets["local_atac"] = redact_path(atac_file)
        assets["local_atac_present"] = Path(atac_file).exists()

    plan: dict[str, Any] = {
        "tissue_id": tissue_id,
        "lane": "construct",
        "g_atac_tag": tag,
        "compare_to_tags": list(COMPARE_TAGS),
        "outputs": {
            "g_atac": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/G_ATAC_v2_{tag}.npz",
            "mantel": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/mantel_vs_locked.json",
            "decomp": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/additive_decomp_row.json",
            "summary": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/construct_summary.json",
        },
        "peerj_decomp_rows_unchanged": 3,
        "peerj_support_rows_unchanged": 13,
        "assets": assets,
        "env": fibro_style_env(tissue_id, atac_file=atac_file),
        "panel_policy": meta.get("panel_policy"),
        "notes": meta.get("notes"),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
    }

    g_path = resolve_g_atac_npz(tag)
    if g_path is None:
        plan["status"] = "awaiting_g_atac"
        plan["next_steps"] = [
            "No local G_ATAC NPZ for this tag — Mantel/decomp cannot run yet.",
            "This wave does NOT treat build_atac_graph_v2.py as a success path.",
            "See docs/reports/download_approval_optional_pilots.md for any future "
            "human-approved data step; then place NPZ under results/v2/extension/construct/{tag}/ "
            "or reuse a locked results/v2/G_ATAC_v2_*.npz before --execute.",
        ]
        return plan

    if not execute:
        plan["status"] = "ready"
        plan["g_atac_source"] = (
            str(g_path.relative_to(ROOT)) if g_path.is_relative_to(ROOT) else str(g_path)
        )
        plan["next_steps"] = [
            f"Re-run with --execute to write Mantel/decomp under results/v2/extension/construct/{tag}/"
        ]
        return plan

    ga, tf_a = _load_consensus(g_path)
    pairs = []
    for other in COMPARE_TAGS:
        if other == tag:
            continue
        other_path = resolve_g_atac_npz(other)
        if other_path is None:
            pairs.append({"pair": [tag, other], "status": "missing_proxy", "path": None})
            continue
        gb, tf_b = _load_consensus(other_path)
        row = _pair_mantel_decomp(tag, ga, tf_a, other, gb, tf_b)
        row["status"] = "ok"
        row["proxy_path"] = str(other_path.relative_to(ROOT))
        pairs.append(row)

    out_dir = EXTENSION_ROOT / "construct" / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    mantel_doc = {
        "schema_version": 1,
        "tissue_id": tissue_id,
        "g_atac_tag": tag,
        "g_atac_source": str(g_path.relative_to(ROOT)) if g_path.is_relative_to(ROOT) else str(g_path),
        "peerj_support_rows_touched": False,
        "pairs": pairs,
    }
    decomp_doc = {
        "schema_version": 1,
        "analysis": "extension_construct_additive_decomp",
        "source_tag": tag,
        "rows": [p for p in pairs if p.get("status") == "ok"],
        "peerj_cross_tissue_additive_decomp_unchanged": True,
    }
    (out_dir / "mantel_vs_locked.json").write_text(json.dumps(mantel_doc, indent=2) + "\n")
    (out_dir / "additive_decomp_row.json").write_text(json.dumps(decomp_doc, indent=2) + "\n")

    plan["status"] = "executed"
    plan["g_atac_source"] = mantel_doc["g_atac_source"]
    plan["n_pairs_ok"] = sum(1 for p in pairs if p.get("status") == "ok")
    plan["mantel_pairs"] = pairs
    (out_dir / "construct_summary.json").write_text(json.dumps(plan, indent=2) + "\n")
    return plan


def mantel_decomp_plan(tissue_id: str) -> dict[str, Any]:
    """Backward-compatible dry-run plan (no execute)."""
    return run_construct(tissue_id, execute=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissue", default="fibroblast", help="Registry tissue id")
    parser.add_argument("--atac-file", default="", help="Optional explicit ATAC/h5ad path")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Compute Mantel/decomp when local G_ATAC exists",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write dry-run / summary JSON under docs/reports/extension-plans/",
    )
    args = parser.parse_args(argv)

    plan = run_construct(
        args.tissue,
        atac_file=args.atac_file or None,
        execute=args.execute,
    )
    print(json.dumps(plan, indent=2))

    if args.write:
        out = PLAN_ROOT / "construct" / str(plan["g_atac_tag"])
        out.mkdir(parents=True, exist_ok=True)
        name = "construct_plan.executed.json" if args.execute else "construct_plan.dry_run.json"
        path = out / name
        path.write_text(json.dumps(plan, indent=2) + "\n")
        print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
