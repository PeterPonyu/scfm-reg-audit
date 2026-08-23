"""Resolve local extension / PeerJ freeze artifact paths (no downloads)."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

EXT_DIR = Path(__file__).resolve().parent
ROOT = EXT_DIR.parents[2]

# PeerJ-safe overlay roots (never write into frozen results/*.public.json).
# results/v2/ is already a LOCAL_WORKTREE_PREFIX in validate_artifacts.py.
EXTENSION_ROOT = ROOT / "results" / "v2" / "extension"
HEAVY_ARTIFACT_ROOT = "results/v2/extension/"
CLAIM_PACK_ROOT = "docs/reports/extension-claim-pack/"
PLAN_ROOT = ROOT / "docs" / "reports" / "extension-plans"

# Locked construct-valid proxies (read-only for Mantel/decomp).
LOCKED_G_ATAC = {
    "GSE174367": ROOT / "results" / "v2" / "G_ATAC_v2_GSE174367.npz",
    "PBMC10k": ROOT / "results" / "v2" / "G_ATAC_v2_PBMC10k.npz",
    "GSE206767": ROOT / "results" / "v2" / "G_ATAC_v2_GSE206767.npz",
}

DESKTOP_DATA = Path(
    os.environ.get("DESKTOP_DATA")
    or os.environ.get("SCREG_DATA_ROOT")
    or (Path.home() / "Desktop" / "data")
)

LOCAL_ATAC_HINTS: dict[str, list[Path]] = {
    "fibroblast": [
        DESKTOP_DATA / "datasets" / "ATAC_data" / "GSE206767_filtered_peak_bc_matrix.h5ad",
        DESKTOP_DATA / "external" / "scfm-reg-audit" / "gse206767",
    ],
    "bmmc": [
        DESKTOP_DATA
        / "external"
        / "scfm-reg-audit"
        / "gse194122"
        / "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad",
    ],
    "brain": [
        DESKTOP_DATA
        / "datasets"
        / "ATAC_data"
        / "GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad",
    ],
}


def resolve_local_atac(tissue_id: str, explicit: str | None = None) -> Path | None:
    """Return first existing local ATAC/h5ad path for a tissue, else None."""
    if explicit:
        p = Path(explicit)
        return p if p.exists() else None
    env_key = {
        "bmmc": "SCREG_BMMC_H5AD",
        "fibroblast": "SCREG_FIBRO_ATAC",
        "brain": "SCREG_BRAIN_ATAC",
    }.get(tissue_id)
    if env_key and os.environ.get(env_key):
        p = Path(os.environ[env_key])
        if p.exists():
            return p
    for p in LOCAL_ATAC_HINTS.get(tissue_id, []):
        if p.exists():
            return p
    return None


def resolve_g_atac_npz(tag: str, prefer_extension: bool = True) -> Path | None:
    """Locate a G_ATAC NPZ by tag (extension overlay first, then locked v2)."""
    candidates: list[Path] = []
    if prefer_extension:
        candidates.append(EXTENSION_ROOT / "construct" / tag / f"G_ATAC_v2_{tag}.npz")
        candidates.append(ROOT / "results" / "v2" / "extension" / "construct" / tag / f"G_ATAC_v2_{tag}.npz")
    locked = LOCKED_G_ATAC.get(tag)
    if locked is not None:
        candidates.append(locked)
    candidates.append(ROOT / "results" / "v2" / f"G_ATAC_v2_{tag}.npz")
    for p in candidates:
        if p.exists():
            return p
    return None


def redact_path(path: Path | str | None) -> str | None:
    """Return repo-relative or ${DESKTOP_DATA}-relative path (no home leaks)."""
    if path is None:
        return None
    p = Path(path)
    try:
        if p.is_relative_to(ROOT):
            return p.relative_to(ROOT).as_posix()
    except (ValueError, OSError):
        pass
    try:
        if p.is_relative_to(DESKTOP_DATA):
            return "${DESKTOP_DATA}/" + p.relative_to(DESKTOP_DATA).as_posix()
    except (ValueError, OSError):
        pass
    # Last resort: basename only (never emit /home/...)
    return p.name


def local_asset_report(tissue_id: str, g_atac_tag: str) -> dict[str, Any]:
    atac = resolve_local_atac(tissue_id)
    g_atac = resolve_g_atac_npz(g_atac_tag)
    return {
        "tissue_id": tissue_id,
        "g_atac_tag": g_atac_tag,
        "local_atac": redact_path(atac),
        "local_atac_present": atac is not None,
        "g_atac_npz": redact_path(g_atac),
        "g_atac_present": g_atac is not None,
        "extension_root": HEAVY_ARTIFACT_ROOT,
        "desktop_data": "${DESKTOP_DATA}",
    }
