#!/usr/bin/env python3
"""D1 DESCARTES spleen RDS → builder-ready ATAC_FILE bridge (no network fetch).

Fail-closed when local RDS/h5ad are absent. Conversion of RDS requires local R
+ anndata tooling and is never attempted via download. Heavy outputs land under
``results/v2/extension/`` only when explicitly prepared.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from paths import (  # noqa: E402
    DESCARTES_PILOT_DIR,
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

TISSUE_ID = "descartes_spleen"
G_ATAC_TAG = "DESCARTES_spleen"
EXPECTED_H5AD_NAME = "descartes_spleen_peaks.h5ad"
RDS_GLOBS = ("*.rds", "*.RDS")


def expected_layout() -> dict[str, str]:
    return {
        "pilot_dir": "${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/",
        "preferred_h5ad": (
            "${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/"
            + EXPECTED_H5AD_NAME
        ),
        "accepted_rds": "${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/*.rds",
        "extension_out": f"{HEAVY_ARTIFACT_ROOT}construct/{G_ATAC_TAG}/",
        "prepared_h5ad_overlay": (
            f"{HEAVY_ARTIFACT_ROOT}construct/{G_ATAC_TAG}/{EXPECTED_H5AD_NAME}"
        ),
    }


def _list_rds(pilot_dir: Path) -> list[Path]:
    if not pilot_dir.is_dir():
        return []
    found: list[Path] = []
    for pat in RDS_GLOBS:
        found.extend(sorted(pilot_dir.glob(pat)))
    # de-dupe while preserving order
    seen: set[Path] = set()
    out: list[Path] = []
    for p in found:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(p)
    return out


def _candidate_h5ads(pilot_dir: Path) -> list[Path]:
    cands: list[Path] = []
    preferred = pilot_dir / EXPECTED_H5AD_NAME
    if preferred.exists():
        cands.append(preferred)
    overlay = EXTENSION_ROOT / "construct" / G_ATAC_TAG / EXPECTED_H5AD_NAME
    if overlay.exists():
        cands.append(overlay)
    if pilot_dir.is_dir():
        for p in sorted(pilot_dir.glob("*.h5ad")):
            if p not in cands:
                cands.append(p)
    return cands


def validate_peak_h5ad(path: Path, *, sample: int = 32) -> dict[str, Any]:
    """Lightweight peak-name check (does not load full matrix into dense)."""
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover
        return {
            "ok": False,
            "path": redact_path(path),
            "error": f"anndata not available: {exc}",
        }
    try:
        A = ad.read_h5ad(path, backed="r")
    except Exception as exc:  # noqa: BLE001 — surface convert/IO errors
        return {
            "ok": False,
            "path": redact_path(path),
            "error": f"failed to open h5ad: {exc}",
        }
    names = [str(x) for x in list(A.var_names[:sample])]
    n_ready = sum(1 for n in names if peak_name_is_builder_ready(n))
    n_fixable = sum(
        1
        for n in names
        if (not peak_name_is_builder_ready(n))
        and peak_name_is_builder_ready(normalize_peak_name(n))
    )
    ok = n_ready == len(names) and len(names) > 0
    return {
        "ok": ok,
        "path": redact_path(path),
        "n_obs": int(A.n_obs),
        "n_vars": int(A.n_vars),
        "sample_var_names": names[:8],
        "sample_builder_ready": n_ready,
        "sample_hyphen_fixable": n_fixable,
        "needs_peak_rename": n_fixable > 0 and n_ready < len(names),
        "note": (
            "Peak names must be chrom:start-end for build_atac_graph_v2.py"
            if not ok
            else "h5ad looks builder-ready on sampled peaks"
        ),
    }


def bridge_status(*, pilot_dir: Path | None = None) -> dict[str, Any]:
    """Report D1 local readiness. Never downloads."""
    pilot = pilot_dir or DESCARTES_PILOT_DIR
    rds = _list_rds(pilot)
    h5ads = _candidate_h5ads(pilot)
    validations = [validate_peak_h5ad(p) for p in h5ads]
    ready = next((v for v in validations if v.get("ok")), None)
    out_rel = f"{HEAVY_ARTIFACT_ROOT}construct/{G_ATAC_TAG}"
    status: dict[str, Any] = {
        "schema_version": 1,
        "plan_id": "D1",
        "tissue_id": TISSUE_ID,
        "g_atac_tag": G_ATAC_TAG,
        "network_fetch_performed": False,
        "pilot_dir": redact_path(pilot) or str(pilot),
        "pilot_dir_present": pilot.is_dir(),
        "rds_files": [redact_path(p) for p in rds],
        "h5ad_candidates": [redact_path(p) for p in h5ads],
        "h5ad_validations": validations,
        "expected_layout": expected_layout(),
        "extension_out": out_rel,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "desktop_data": "${DESKTOP_DATA}",
        "desktop_data_resolved_redacted": redact_path(DESKTOP_DATA),
    }
    if ready is not None:
        atac = next(p for p, v in zip(h5ads, validations) if v is ready)
        status["status"] = "ready_atac_file"
        status["ATAC_FILE"] = str(atac.resolve())
        status["ATAC_FILE_REDACTED"] = redact_path(atac)
        status["build_command"] = builder_env_command(
            tag=G_ATAC_TAG,
            atac_file=str(atac.resolve()),
            extension_out=out_rel,
            meta_file="none",
        )
        status["next_steps"] = [
            "Local builder-ready h5ad present — no download needed.",
            f"Build overlay G_ATAC: {status['build_command']}",
            (
                "Then: python src/v2/extension/cli.py construct "
                f"--tissue {TISSUE_ID} --execute --write"
            ),
        ]
        return status

    if rds:
        status["status"] = "rds_present_needs_conversion"
        status["message"] = (
            f"Found RDS under {redact_path(pilot)} but no builder-ready h5ad. "
            "Convert locally (R/Seurat or Matrix → peak matrix h5ad with "
            "var_names chrom:start-end); this CLI does not fetch or run R."
        )
        status["next_steps"] = [
            f"Place converted peaks h5ad at {expected_layout()['preferred_h5ad']}",
            "Or write under " + expected_layout()["prepared_h5ad_overlay"],
            "Re-run: python src/v2/extension/cli.py descartes-bridge",
            "Then build with SCREG_EXTENSION_OUT + SCREG_PEERJ_SUPPORT_LOCK=1",
        ]
        status["r_conversion_sketch"] = [
            "# Outside this CLI (no network). Example sketch only:",
            "# R: readRDS → extract peak-by-cell counts + peak coords",
            "# write Matrix Market / CSV, then python anndata write_h5ad",
            f"# Target: {expected_layout()['preferred_h5ad']}",
        ]
        return status

    status["status"] = "absent_local"
    status["message"] = (
        "No DESCARTES spleen RDS/h5ad on this machine. "
        f"After human-approved D1 fetch, place files under "
        f"{expected_layout()['pilot_dir']} "
        f"(preferred h5ad: {expected_layout()['preferred_h5ad']}). "
        "This bridge never downloads."
    )
    status["next_steps"] = [
        "Fill docs/reports/download_approval_optional_pilots.md (D1) + env gates",
        "Manual fetch only after approval (not via extension CLI)",
        f"Place RDS/h5ad under {expected_layout()['pilot_dir']}",
        "Re-run descartes-bridge → build_atac_graph_v2 with SCREG_EXTENSION_OUT",
    ]
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot-dir",
        default="",
        help="Override local DESCARTES pilot directory (no download)",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write status JSON under results/v2/extension/construct/DESCARTES_spleen/",
    )
    args = parser.parse_args(argv)
    pilot = Path(args.pilot_dir) if args.pilot_dir else None
    status = bridge_status(pilot_dir=pilot)
    print(json.dumps(status, indent=2))
    if args.write:
        out_dir = assert_confined_write_path(
            EXTENSION_ROOT / "construct" / G_ATAC_TAG,
            label="descartes bridge out",
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "descartes_bridge_status.json"
        # PeerJ-safe: drop absolute ATAC_FILE from written artifact
        doc = dict(status)
        if "ATAC_FILE" in doc:
            doc["ATAC_FILE"] = doc.get("ATAC_FILE_REDACTED") or doc["ATAC_FILE"]
        path.write_text(json.dumps(doc, indent=2) + "\n")
        print(f"wrote {path.relative_to(ROOT)}", file=sys.stderr)
    # Exit 0 when ready; 2 when absent / needs conversion (fail-closed readiness).
    if status.get("status") == "ready_atac_file":
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
