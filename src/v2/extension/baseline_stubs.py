#!/usr/bin/env python3
"""Runnable baseline emitters for extension Tier A–C comparisons.

Uses local locked ``G_ATAC`` / public JSON only. No FM BH membership writes.
No PeerJ Support mutation. Outputs under ``results/v2/extension/baselines/``.
"""

from __future__ import annotations

import argparse
import json
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
    redact_path,
    resolve_g_atac_npz,
)
from registry import load_extension_registry  # noqa: E402

STUB_METHODS = (
    "degree_matched_random",
    "motif_only_rp",
    "encode_chip_binding",
    "collectri_prior",
)

DEFAULT_PROXY_TAG = "GSE174367"


def _load_consensus(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    z = np.load(path, allow_pickle=False)
    types = [str(t) for t in z["types"]]
    g = np.mean([z[f"G_{t}"] for t in types], axis=0).astype(np.float64)
    genes = [str(gname) for gname in z["genes"]]
    return g, np.asarray(z["tf_rows"], dtype=int), genes


def _tf_submatrix(g: np.ndarray, tf_rows: np.ndarray) -> np.ndarray:
    sub = g[tf_rows].copy()
    for k, tf in enumerate(tf_rows):
        sub[k, int(tf)] = 0.0
    return sub


def emit_degree_matched_random(proxy_tag: str = DEFAULT_PROXY_TAG, seed: int = 0) -> dict[str, Any]:
    path = resolve_g_atac_npz(proxy_tag)
    if path is None:
        return {"method_id": "degree_matched_random", "status": "skipped_missing_g_atac", "proxy_tag": proxy_tag}
    g, tf_rows, _genes = _load_consensus(path)
    sub = _tf_submatrix(g, tf_rows)
    out_deg = (sub > 0).sum(axis=1)
    in_deg = (sub > 0).sum(axis=0)
    rng = np.random.default_rng(seed)
    # Configuration-model style: sample random TF→gene edges matching out-degree,
    # then score with Uniform(0,1) weights (topology null, not FM).
    scores = np.zeros_like(sub)
    n_genes = sub.shape[1]
    for i, deg in enumerate(out_deg):
        deg_i = int(deg)
        if deg_i <= 0:
            continue
        # avoid self index if TF gene index is in range
        banned = {int(tf_rows[i])}
        choices = [j for j in range(n_genes) if j not in banned]
        if deg_i > len(choices):
            deg_i = len(choices)
        picks = rng.choice(choices, size=deg_i, replace=False)
        scores[i, picks] = rng.random(deg_i)
    # Agreement of null binary support vs true binary support (should be ~chance).
    true_bin = (sub > 0).astype(np.float64).ravel()
    null_bin = (scores > 0).astype(np.float64).ravel()
    if true_bin.std() > 0 and null_bin.std() > 0:
        phi = float(np.corrcoef(true_bin, null_bin)[0, 1])
    else:
        phi = None
    return {
        "method_id": "degree_matched_random",
        "status": "emitted",
        "proxy_tag": proxy_tag,
        "proxy_path": str(path.relative_to(ROOT)),
        "seed": seed,
        "shape": [int(scores.shape[0]), int(scores.shape[1])],
        "n_edges": int((scores > 0).sum()),
        "mean_out_degree": float(out_deg.mean()),
        "mean_in_degree_targets": float(in_deg.mean()),
        "binary_phi_vs_proxy": None if phi is None else round(phi, 6),
        "score_matrix_relpath": f"{HEAVY_ARTIFACT_ROOT}baselines/degree_matched_random/scores.npz",
        "scores": scores,
        "tf_rows": tf_rows,
    }


def emit_motif_only_rp(proxy_tag: str = DEFAULT_PROXY_TAG) -> dict[str, Any]:
    """Ablate ATAC magnitude: binary motif/support mask from G_ATAC > 0."""
    path = resolve_g_atac_npz(proxy_tag)
    if path is None:
        return {"method_id": "motif_only_rp", "status": "skipped_missing_g_atac", "proxy_tag": proxy_tag}
    g, tf_rows, _genes = _load_consensus(path)
    sub = _tf_submatrix(g, tf_rows)
    binary = (sub > 0).astype(np.float64)
    # Compare weighted scores vs binary mask over full TF×gene submatrix
    # (on positive edges alone binary is constant → Spearman undefined).
    if sub.size < 10 or float(np.std(sub)) == 0.0 or float(np.std(binary)) == 0.0:
        rho = None
    else:
        rho = float(spearmanr(sub.ravel(), binary.ravel()).statistic)
    return {
        "method_id": "motif_only_rp",
        "status": "emitted",
        "proxy_tag": proxy_tag,
        "proxy_path": str(path.relative_to(ROOT)),
        "shape": [int(binary.shape[0]), int(binary.shape[1])],
        "n_edges": int(binary.sum()),
        "spearman_weighted_vs_binary_on_edges": None if rho is None else round(rho, 6),
        "note": "Binary support from local G_ATAC; ATAC magnitude ablated (simple-method SI).",
        "score_matrix_relpath": f"{HEAVY_ARTIFACT_ROOT}baselines/motif_only_rp/scores.npz",
        "scores": binary,
        "tf_rows": tf_rows,
    }


def emit_encode_chip_binding() -> dict[str, Any]:
    pub = ROOT / "results" / "encode_proxy_calibration_v1.public.json"
    meta = (
        Path.home()
        / "Desktop"
        / "data"
        / "datasets"
        / "extension_pilots"
        / "encode_chip"
        / "encode_tf_chip_human_hg38_metadata.json"
    )
    if not pub.exists():
        return {"method_id": "encode_chip_binding", "status": "skipped_missing_public_json"}
    doc = json.loads(pub.read_text())
    tissues = doc.get("tissues") or {}
    summary_rows = []
    if isinstance(tissues, dict):
        for tid, payload in tissues.items():
            if isinstance(payload, dict):
                summary_rows.append(
                    {
                        "tissue": tid,
                        "keys": sorted(payload.keys())[:12],
                        "n_keys": len(payload),
                    }
                )
    return {
        "method_id": "encode_chip_binding",
        "status": "emitted",
        "source_public_json": "results/encode_proxy_calibration_v1.public.json",
        "local_metadata_present": meta.exists(),
        "local_metadata": redact_path(meta) if meta.exists() else None,
        "schema_version": doc.get("schema_version"),
        "job": doc.get("job"),
        "coverage_gate": doc.get("coverage_gate"),
        "n_panel_genes": doc.get("n_panel_genes"),
        "n_panel_tfs": doc.get("n_panel_tfs"),
        "tissue_summaries": summary_rows,
        "note": "Reuses shipped ENCODE calibration JSON; no new download.",
    }


def emit_collectri_prior() -> dict[str, Any]:
    candidates = [
        Path.home() / "Desktop" / "data" / "datasets" / "extension_pilots" / "collectri",
        ROOT / "data" / "collectri",
        ROOT / "results" / "v2" / "extension" / "baselines" / "collectri_prior" / "cache",
    ]
    found = next((p for p in candidates if p.exists()), None)
    if found is None:
        return {
            "method_id": "collectri_prior",
            "status": "skipped_no_local_cache",
            "searched": [redact_path(p) for p in candidates],
            "note": "No CollecTRI/OmniPath cache locally — approval-only until cache present.",
        }
    files = sorted(found.rglob("*"))[:20]
    return {
        "method_id": "collectri_prior",
        "status": "cache_present_not_parsed",
        "cache_root": redact_path(found),
        "n_files_seen": len(files),
        "sample_files": [
            f.relative_to(found).as_posix() if f.is_relative_to(found) else f.name
            for f in files[:8]
        ],
        "note": "Cache detected; full prior→panel projection left for a follow-up parse job.",
    }


def baseline_plan(method_id: str) -> dict[str, Any]:
    """Dry-run metadata plan (no score matrices)."""
    reg = load_extension_registry()
    meta = reg.get_method(method_id)
    if meta.get("bh_membership") == "fm_families":
        raise PermissionError(f"{method_id} is an FM-family method, not a baseline stub")
    return {
        "method_id": method_id,
        "tier": meta.get("tier"),
        "estimand": meta.get("estimand"),
        "bh_membership": meta.get("bh_membership"),
        "status": meta.get("status", "stub"),
        "panel": "frozen_446x1200",
        "outputs": {
            "score_matrix": f"{HEAVY_ARTIFACT_ROOT}baselines/{method_id}/scores.npz",
            "summary": f"{HEAVY_ARTIFACT_ROOT}baselines/{method_id}/summary.json",
        },
        "peerj_bh_families_unchanged": 8,
        "hooks": {
            "degree_matched_random": "degree-matched random edges on panel from local G_ATAC",
            "motif_only_rp": "binary G_ATAC support (ATAC magnitude ablated)",
            "encode_chip_binding": "summarize results/encode_proxy_calibration_v1.public.json",
            "collectri_prior": "use local CollecTRI cache if present; else skip",
        }.get(method_id, "implement score matrix on panel"),
        "notes": meta.get("notes"),
    }


def run_baseline(method_id: str, *, execute: bool = False, proxy_tag: str = DEFAULT_PROXY_TAG) -> dict[str, Any]:
    plan = baseline_plan(method_id)
    if not execute:
        plan["run_status"] = "dry_run"
        return plan

    emitters = {
        "degree_matched_random": lambda: emit_degree_matched_random(proxy_tag=proxy_tag),
        "motif_only_rp": lambda: emit_motif_only_rp(proxy_tag=proxy_tag),
        "encode_chip_binding": emit_encode_chip_binding,
        "collectri_prior": emit_collectri_prior,
    }
    payload = emitters[method_id]()
    out_dir = EXTENSION_ROOT / "baselines" / method_id
    out_dir.mkdir(parents=True, exist_ok=True)
    scores = payload.pop("scores", None)
    tf_rows = payload.pop("tf_rows", None)
    if scores is not None:
        np.savez_compressed(
            out_dir / "scores.npz",
            scores=scores.astype(np.float32),
            tf_rows=np.asarray(tf_rows, dtype=np.int32),
            proxy_tag=np.asarray(proxy_tag),
        )
    summary = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "peerj_support_rows_touched": False,
        "peerj_bh_families_unchanged": 8,
        "plan": plan,
        "result": payload,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plan["run_status"] = payload.get("status", "emitted")
    plan["result"] = payload
    plan["summary_path"] = str((out_dir / "summary.json").relative_to(ROOT))
    return plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=STUB_METHODS, default="motif_only_rp")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--execute", action="store_true", help="Emit real baseline artifacts")
    parser.add_argument("--proxy-tag", default=DEFAULT_PROXY_TAG)
    parser.add_argument("--write", action="store_true", help="Also write plan JSON under docs/reports")
    args = parser.parse_args(argv)

    methods = list(STUB_METHODS) if args.all else [args.method]
    plans = [
        run_baseline(mid, execute=args.execute, proxy_tag=args.proxy_tag) for mid in methods
    ]
    printable = []
    for p in plans:
        # Avoid dumping huge nested plan twice
        printable.append({k: v for k, v in p.items() if k != "result" or args.execute})
    print(json.dumps(printable if args.all else printable[0], indent=2))

    if args.write:
        for plan in plans:
            out = PLAN_ROOT / "baselines" / plan["method_id"]
            out.mkdir(parents=True, exist_ok=True)
            name = (
                "baseline_plan.executed.json" if args.execute else "baseline_plan.dry_run.json"
            )
            slim = {k: v for k, v in plan.items() if k != "result"}
            if "result" in plan:
                slim["result"] = plan["result"]
            path = out / name
            path.write_text(json.dumps(slim, indent=2) + "\n")
            print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
