#!/usr/bin/env python3
"""D3 HTAN GBM single-sample pilot: local tar inventory → prepare or blocked.

Does **not** download. Does **not** call peaks (MACS2/Signac). Default path for
the on-disk Cell Ranger fragments-only tar is structured ``status=blocked``
(wave dual-path E4.3b). Happy path only if a builder-ready peak h5ad already
exists under the pilot dir or extension overlay.
"""

from __future__ import annotations

import argparse
import json
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from paths import (  # noqa: E402
    DESKTOP_DATA,
    EXTENSION_ROOT,
    HEAVY_ARTIFACT_ROOT,
    ROOT,
    assert_confined_write_path,
    builder_env_command,
    normalize_peak_name,
    peak_name_is_builder_ready,
    redact_path,
)

TISSUE_ID = "htan_gbm_pilot"
G_ATAC_TAG = "HTAN_GBM_C3N01334"
PLAN_ID = "D3"
DEFAULT_PILOT_DIR = (
    DESKTOP_DATA / "datasets" / "extension_pilots" / "htan" / "sample_pilot"
)
DEFAULT_TAR_NAME = "GSM7710026_C3N-01334_CPT0125220004_snATAC_GBM.tar.gz"
PREPARED_ATAC_NAME = f"{G_ATAC_TAG}_atac_peaks.h5ad"
PREPARED_META_NAME = f"{G_ATAC_TAG}_cell_meta.csv.gz"
STATUS_NAME = "htan_prepare_status.json"


def prepared_paths(tag: str = G_ATAC_TAG) -> dict[str, Path]:
    out = EXTENSION_ROOT / "construct" / tag
    return {
        "out_dir": out,
        "atac": out / PREPARED_ATAC_NAME,
        "meta": out / PREPARED_META_NAME,
        "status": out / STATUS_NAME,
    }


def resolve_tar(pilot_dir: Path, tar: Path | None = None) -> Path | None:
    if tar is not None:
        return tar if tar.exists() else None
    exact = pilot_dir / DEFAULT_TAR_NAME
    if exact.exists():
        return exact
    matches = sorted(pilot_dir.glob("GSM7710026_*_snATAC_GBM.tar.gz"))
    if matches:
        return matches[0]
    matches = sorted(pilot_dir.glob("*.tar.gz"))
    return matches[0] if matches else None


def inventory_htan_tar(tar_path: Path) -> dict[str, Any]:
    """Classify tar members without full extract (tarfile listing only)."""
    members: list[str] = []
    with tarfile.open(tar_path, "r:gz") as tf:
        for m in tf.getmembers():
            members.append(m.name)

    lower = [m.lower() for m in members]
    has_fragments = any("fragment" in m and m.endswith(".gz") for m in lower) or any(
        m.endswith("fragments.tsv.gz") for m in lower
    )
    has_peaks = any(
        ("peak" in m and (m.endswith(".bed") or m.endswith(".bed.gz") or "peak_annotation" in m))
        for m in lower
    )
    has_matrix = any(
        any(tok in m for tok in ("peak_bc", "filtered_peak", ".h5", "matrix.mtx", "barcodes.tsv"))
        for m in lower
    )
    # pure path tokens for h5/mtx are weak; refine
    has_peak_bc_matrix = any(
        "peak" in m and any(x in m for x in (".h5", "mtx", "matrix"))
        for m in lower
    ) or any("filtered_peak_bc_matrix" in m for m in lower)
    has_cell_type_meta = any(
        any(x in m for x in ("cell_type", "celltype", "singlecell.csv", "meta"))
        and not m.endswith(".tar.gz")
        for m in lower
    )
    fragments_member = next(
        (m for m in members if m.lower().endswith("fragments.tsv.gz")),
        None,
    )

    if has_fragments and not has_peaks and not has_peak_bc_matrix and not has_matrix:
        layout = "fragments_only"
    elif has_peak_bc_matrix or has_peaks:
        layout = "peak_matrix_or_peaks_present"
    else:
        layout = "unknown"

    return {
        "n_members": len(members),
        "members_sample": members[:20],
        "has_fragments": has_fragments,
        "has_peaks": has_peaks,
        "has_peak_bc_matrix": has_peak_bc_matrix or has_matrix,
        "has_cell_type_meta": has_cell_type_meta,
        "fragments_member": fragments_member,
        "layout": layout,
        "pipeline": "cellranger-atac-2.0.0" if has_fragments else None,
        "genome": "hg38",
    }


def _find_local_peak_h5ad(pilot_dir: Path, paths: dict[str, Path]) -> Path | None:
    candidates = [
        paths["atac"],
        pilot_dir / PREPARED_ATAC_NAME,
        pilot_dir / "atac_peaks.h5ad",
        pilot_dir / "unpacked" / PREPARED_ATAC_NAME,
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def _inspect_peak_h5ad(path: Path, *, sample: int = 8) -> dict[str, Any]:
    import anndata as ad

    A = ad.read_h5ad(path, backed="r")
    names = [str(x) for x in list(A.var_names[:sample])]
    ready = sum(1 for n in names if peak_name_is_builder_ready(n))
    renamable = sum(
        1
        for n in names
        if (not peak_name_is_builder_ready(n))
        and peak_name_is_builder_ready(normalize_peak_name(n))
    )
    return {
        "path": redact_path(path),
        "n_obs": int(A.n_obs),
        "n_vars": int(A.n_vars),
        "sample_var_names": names,
        "sample_builder_ready": ready,
        "sample_needs_rename": renamable,
        "builder_ready": ready == len(names) and len(names) > 0,
    }


def prepare_status(
    *,
    pilot_dir: Path | None = None,
    tar: Path | None = None,
    execute: bool = False,
) -> dict[str, Any]:
    """Dry-run or execute prepare; fragments-only tar → structured blocked."""
    del execute  # no peak-calling; execute does not invent matrices
    pilot = pilot_dir or DEFAULT_PILOT_DIR
    paths = prepared_paths()
    out_rel = f"{HEAVY_ARTIFACT_ROOT}construct/{G_ATAC_TAG}"
    assert_confined_write_path(out_rel, label="SCREG_EXTENSION_OUT")

    doc: dict[str, Any] = {
        "schema_version": 1,
        "plan_id": PLAN_ID,
        "tissue_id": TISSUE_ID,
        "g_atac_tag": G_ATAC_TAG,
        "role": "construct_candidate",
        "network_fetch_performed": False,
        "fm_support_allowed": False,
        "panel_policy_gate": "construct_candidate_only",
        "extension_out": out_rel,
        "pilot_dir": redact_path(pilot),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "peak_set_strategies_considered": [
            "call_peaks",
            "external_peaks",
            "fixed_bins",
            "blocked",
        ],
    }

    tar_path = resolve_tar(pilot, tar)
    doc["tar_present"] = tar_path is not None and tar_path.exists()
    if tar_path is not None:
        doc["tar_path"] = redact_path(tar_path)
        try:
            doc["tar_bytes"] = int(tar_path.stat().st_size)
        except OSError:
            doc["tar_bytes"] = None
    else:
        doc["tar_path"] = None
        doc["tar_bytes"] = None

    # Happy path: pre-existing builder-ready peak matrix on disk.
    peak_h5ad = _find_local_peak_h5ad(pilot, paths)
    if peak_h5ad is not None:
        try:
            insp = _inspect_peak_h5ad(peak_h5ad)
        except Exception as exc:  # noqa: BLE001
            doc["status"] = "inspect_failed"
            doc["message"] = str(exc)
            return doc
        doc["inspect"] = insp
        if insp.get("builder_ready") or insp.get("sample_needs_rename", 0) == len(
            insp.get("sample_var_names") or []
        ):
            meta = paths["meta"] if paths["meta"].exists() else None
            doc["status"] = "prepared_present" if peak_h5ad == paths["atac"] else "ready_to_prepare"
            doc["peak_set_strategy"] = "external_or_prebuilt_h5ad"
            doc["prepared_atac"] = redact_path(peak_h5ad)
            doc["build_command"] = builder_env_command(
                tag=G_ATAC_TAG,
                atac_file=str(peak_h5ad.resolve()),
                extension_out=out_rel,
                meta_file=str(meta.resolve()) if meta else "none",
            )
            doc["accepted_wave_completion"] = "E4.3a"
            doc["next_steps"] = [
                f"Build overlay G_ATAC via emitted build_command (lock+OUT).",
                f"Then: python src/v2/extension/cli.py construct --tissue {TISSUE_ID} --execute --write",
                "Never inflate PeerJ 13-row SAP / FM Support.",
            ]
            return doc

    if not doc["tar_present"]:
        doc["status"] = "absent_local_tar"
        doc["message"] = (
            "No HTAN pilot tar under pilot_dir. Place GSM7710026_*_snATAC_GBM.tar.gz "
            "under ${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/ "
            "(no download in this CLI)."
        )
        doc["peak_set_strategy"] = "blocked"
        doc["build_command"] = None
        doc["next_steps"] = [
            "Do not use download --plan-id D3 for fetch; local unpack/prepare only.",
            "When tar is present, re-run: python src/v2/extension/cli.py prepare-htan --write",
        ]
        return doc

    assert tar_path is not None
    try:
        inv = inventory_htan_tar(tar_path)
    except Exception as exc:  # noqa: BLE001
        doc["status"] = "inspect_failed"
        doc["message"] = str(exc)
        return doc
    doc["inventory"] = inv

    if inv.get("layout") == "fragments_only" or (
        inv.get("has_fragments")
        and not inv.get("has_peaks")
        and not inv.get("has_peak_bc_matrix")
    ):
        doc["status"] = "blocked"
        doc["block_reason"] = "fragments_only_no_peak_matrix"
        doc["peak_set_strategy"] = "blocked"
        doc["accepted_wave_completion"] = "E4.3b"
        doc["build_command"] = None
        doc["message"] = (
            "Local D3 tar is Cell Ranger fragments-only; no peaks/matrix/cell-type "
            "meta in archive. Blocker path is intentional wave completion (E4.3b)."
        )
        doc["next_steps"] = [
            "Do not call peaks in default prepare-htan.",
            "Optional future: supply external peak h5ad with chrom:start-end var_names "
            f"as {PREPARED_ATAC_NAME} under pilot_dir or extension overlay.",
            "D3 download gate: decision=no_new_download (asset already local).",
            "Never inflate PeerJ 13-row SAP / FM Support for this tissue.",
        ]
        return doc

    # Non-fragments layout without usable h5ad — still blocked honestly.
    doc["status"] = "blocked"
    doc["block_reason"] = "no_builder_ready_peak_matrix"
    doc["peak_set_strategy"] = "blocked"
    doc["accepted_wave_completion"] = "E4.3b"
    doc["build_command"] = None
    doc["message"] = (
        f"Tar layout={inv.get('layout')} but no builder-ready peak h5ad on disk."
    )
    doc["next_steps"] = [
        "Provide peak matrix h5ad with chrom:start-end var_names, then re-run prepare-htan.",
        "Do not network-fetch D4/D5 in this wave.",
    ]
    return doc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot-dir",
        default="",
        help="Pilot root (default: ${DESKTOP_DATA}/.../htan/sample_pilot/)",
    )
    parser.add_argument("--tar", default="", help="Explicit tar.gz path")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Reserved for happy-path extract only; never auto-calls peaks",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write status JSON under results/v2/extension/construct/HTAN_GBM_C3N01334/",
    )
    args = parser.parse_args(argv)
    pilot = Path(args.pilot_dir) if args.pilot_dir else None
    tar = Path(args.tar) if args.tar else None
    doc = prepare_status(pilot_dir=pilot, tar=tar, execute=args.execute)
    print(json.dumps(doc, indent=2))
    if args.write:
        out_dir = assert_confined_write_path(
            prepared_paths()["out_dir"], label="htan prepare out"
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / STATUS_NAME
        path.write_text(json.dumps(doc, indent=2) + "\n")
        print(f"wrote {path.relative_to(ROOT)}", file=sys.stderr)
    if doc.get("status") in {"inspect_failed"}:
        return 2
    # absent_local_tar is soft (exit 0) when used as inventory; keep 0 for blocked too
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
