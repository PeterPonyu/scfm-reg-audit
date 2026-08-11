#!/usr/bin/env python3
"""Prepare BMMC OpenProblems multiome h5ad → builder-ready ATAC peak matrix.

Does **not** download. Does **not** run full ``build_atac_graph_v2`` (heavy).
Writes peak h5ad + Barcode,Cell.Type meta under ``results/v2/extension/`` when
``--execute`` is set. Default is dry-run / status only.

Peak names ``chr1-start-end`` are rewritten to ``chr1:start-end``.
FM Support / PeerJ 13-row SAP remain out of scope (P3 gate).
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from paths import (  # noqa: E402
    EXTENSION_ROOT,
    HEAVY_ARTIFACT_ROOT,
    ROOT,
    assert_confined_write_path,
    builder_env_command,
    normalize_peak_name,
    peak_name_is_builder_ready,
    redact_path,
    resolve_local_atac,
)

TISSUE_ID = "bmmc"
G_ATAC_TAG = "GSE194122"
PREPARED_ATAC_NAME = "GSE194122_atac_peaks.h5ad"
PREPARED_META_NAME = "GSE194122_cell_meta.csv.gz"


def prepared_paths(tag: str = G_ATAC_TAG) -> dict[str, Path]:
    out = EXTENSION_ROOT / "construct" / tag
    return {
        "out_dir": out,
        "atac": out / PREPARED_ATAC_NAME,
        "meta": out / PREPARED_META_NAME,
    }


def inspect_multiome_h5ad(path: Path, *, sample: int = 16) -> dict[str, Any]:
    """Classify local BMMC/multiome h5ad without building G_ATAC."""
    import anndata as ad

    A = ad.read_h5ad(path, backed="r")
    report: dict[str, Any] = {
        "path": redact_path(path),
        "n_obs": int(A.n_obs),
        "n_vars": int(A.n_vars),
        "has_cell_type": "cell_type" in A.obs.columns,
        "has_feature_types": "feature_types" in A.var.columns,
    }
    if "feature_types" in A.var.columns:
        ft = A.var["feature_types"].astype(str)
        vc = ft.value_counts().to_dict()
        report["feature_type_counts"] = {str(k): int(v) for k, v in vc.items()}
        atac_mask = ft.str.upper().isin(("ATAC", "PEAKS", "PEAK"))
        n_atac = int(atac_mask.sum())
        report["n_atac_features"] = n_atac
        if n_atac > 0:
            names = [str(x) for x in A.var_names[np.where(atac_mask.to_numpy())[0][:sample]]]
        else:
            names = [str(x) for x in list(A.var_names[:sample])]
    else:
        names = [str(x) for x in list(A.var_names[:sample])]
        report["n_atac_features"] = None

    report["sample_var_names"] = names
    report["sample_builder_ready"] = sum(1 for n in names if peak_name_is_builder_ready(n))
    report["sample_needs_rename"] = sum(
        1
        for n in names
        if (not peak_name_is_builder_ready(n))
        and peak_name_is_builder_ready(normalize_peak_name(n))
    )
    # Already a pure peak matrix?
    if report.get("n_atac_features") in (None, report["n_vars"]) and report[
        "sample_builder_ready"
    ] == len(names):
        report["layout"] = "builder_ready_peak_matrix"
    elif report.get("n_atac_features") and report["sample_needs_rename"] > 0:
        report["layout"] = "multiome_atac_needs_extract_and_rename"
    elif report.get("n_atac_features"):
        report["layout"] = "multiome_atac_needs_extract"
    else:
        report["layout"] = "unknown_or_gex_only"
    return report


def extract_atac_peak_matrix(
    src: Path,
    *,
    out_atac: Path,
    out_meta: Path | None = None,
    cell_type_col: str = "cell_type",
) -> dict[str, Any]:
    """Write ATAC-only h5ad with chrom:start-end peaks (+ optional meta CSV)."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    A = ad.read_h5ad(src)
    if "feature_types" in A.var.columns:
        ft = A.var["feature_types"].astype(str).str.upper()
        mask = ft.isin(("ATAC", "PEAKS", "PEAK")).to_numpy()
        if not mask.any():
            raise ValueError(f"no ATAC features in {redact_path(src)}")
        B = A[:, mask].copy()
    else:
        B = A.copy()

    new_names = [normalize_peak_name(str(n)) for n in B.var_names]
    if not all(peak_name_is_builder_ready(n) for n in new_names):
        bad = [n for n in new_names if not peak_name_is_builder_ready(n)][:5]
        raise ValueError(f"peak names not builder-ready after normalize: {bad}")
    B.var_names = pd.Index(new_names)
    B.var_names_make_unique()
    if not sp.issparse(B.X):
        B.X = sp.csr_matrix(B.X)

    out_atac.parent.mkdir(parents=True, exist_ok=True)
    B.write_h5ad(out_atac)

    meta_path = None
    if out_meta is not None and cell_type_col in B.obs.columns:
        out_meta.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(out_meta, "wt") as fh:
            fh.write("Barcode,Cell.Type\n")
            for bc, ct in zip(B.obs_names, B.obs[cell_type_col].astype(str)):
                fh.write(f"{bc},{ct}\n")
        meta_path = out_meta

    return {
        "atac": redact_path(out_atac),
        "meta": redact_path(meta_path) if meta_path else None,
        "n_obs": int(B.n_obs),
        "n_peaks": int(B.n_vars),
    }


def prepare_status(
    *,
    src: Path | None = None,
    execute: bool = False,
) -> dict[str, Any]:
    paths = prepared_paths()
    out_rel = f"{HEAVY_ARTIFACT_ROOT}construct/{G_ATAC_TAG}"
    assert_confined_write_path(out_rel, label="SCREG_EXTENSION_OUT")

    resolved_src = src or resolve_local_atac(TISSUE_ID)
    doc: dict[str, Any] = {
        "schema_version": 1,
        "tissue_id": TISSUE_ID,
        "g_atac_tag": G_ATAC_TAG,
        "panel_policy_gate": "P3_defer_fm_support",
        "fm_support_allowed": False,
        "network_fetch_performed": False,
        "extension_out": out_rel,
        "prepared_atac": redact_path(paths["atac"]),
        "prepared_meta": redact_path(paths["meta"]),
        "prepared_atac_present": paths["atac"].exists(),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
    }

    if resolved_src is None:
        doc["status"] = "absent_local_h5ad"
        doc["message"] = (
            "No local BMMC h5ad. Set SCREG_BMMC_H5AD or place file under "
            "${DESKTOP_DATA}/external/scfm-reg-audit/gse194122/."
        )
        doc["next_steps"] = [
            "No download in this CLI — use already-local D2 asset only.",
            "See docs/reports/bmmc-panel-policy-memo.md (P3: no FM Support).",
        ]
        return doc

    doc["source_h5ad"] = redact_path(resolved_src)
    doc["source_h5ad_present"] = True
    try:
        doc["inspect"] = inspect_multiome_h5ad(resolved_src)
    except Exception as exc:  # noqa: BLE001
        doc["status"] = "inspect_failed"
        doc["message"] = str(exc)
        return doc

    # Prefer already-prepared overlay ATAC if present.
    atac_for_build = paths["atac"] if paths["atac"].exists() else None
    meta_for_build = paths["meta"] if paths["meta"].exists() else None
    if atac_for_build is None and doc["inspect"].get("layout") == "builder_ready_peak_matrix":
        atac_for_build = resolved_src

    if execute and atac_for_build is None:
        try:
            written = extract_atac_peak_matrix(
                resolved_src,
                out_atac=paths["atac"],
                out_meta=paths["meta"],
            )
            doc["extract"] = written
            atac_for_build = paths["atac"]
            meta_for_build = paths["meta"] if paths["meta"].exists() else None
            doc["status"] = "prepared"
        except Exception as exc:  # noqa: BLE001
            doc["status"] = "prepare_failed"
            doc["message"] = str(exc)
            return doc
    elif atac_for_build is not None and paths["atac"].exists():
        doc["status"] = "prepared_present" if not execute else "prepared"
    elif not execute:
        doc["status"] = "ready_to_prepare"
    else:
        doc["status"] = "ready_to_prepare"

    if atac_for_build is not None and atac_for_build.exists():
        meta_arg = str(meta_for_build.resolve()) if meta_for_build and meta_for_build.exists() else "none"
        doc["build_command"] = builder_env_command(
            tag=G_ATAC_TAG,
            atac_file=str(atac_for_build.resolve()),
            extension_out=out_rel,
            meta_file=meta_arg,
        )
        doc["next_steps"] = [
            "P3: construct-lane only — do not inflate PeerJ FM Support rows.",
            f"Build extension G_ATAC (heavy; not run by this prepare step): {doc['build_command']}",
            (
                "Then Mantel/decomp: python src/v2/extension/cli.py construct "
                f"--tissue {TISSUE_ID} --execute --write"
            ),
        ]
    else:
        doc["prepare_command"] = (
            f"python src/v2/extension/cli.py prepare-bmmc --execute"
        )
        doc["next_steps"] = [
            "Local multiome h5ad present but peak ATAC not yet extracted.",
            f"Run: {doc['prepare_command']}",
            "Then run emitted build_command (SCREG_EXTENSION_OUT + PeerJ lock).",
            "P3: no FM Support / no PeerJ 13-row rewrite.",
        ]
    return doc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", default="", help="Explicit multiome h5ad path")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Extract ATAC peaks + meta into results/v2/extension/construct/GSE194122/",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write status JSON next to prepared artifacts",
    )
    args = parser.parse_args(argv)
    src = Path(args.src) if args.src else None
    doc = prepare_status(src=src, execute=args.execute)
    print(json.dumps(doc, indent=2))
    if args.write:
        out_dir = assert_confined_write_path(
            prepared_paths()["out_dir"], label="bmmc prepare out"
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "bmmc_prepare_status.json"
        path.write_text(json.dumps(doc, indent=2) + "\n")
        print(f"wrote {path.relative_to(ROOT)}", file=sys.stderr)
    if doc.get("status") in {
        "absent_local_h5ad",
        "inspect_failed",
        "prepare_failed",
    }:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
