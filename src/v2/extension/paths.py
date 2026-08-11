"""Resolve local extension / PeerJ freeze artifact paths (no downloads)."""

from __future__ import annotations

import os
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

EXT_DIR = Path(__file__).resolve().parent
ROOT = EXT_DIR.parents[2]

# PeerJ-safe overlay roots (never write into frozen results/*.public.json).
# results/v2/ is already a LOCAL_WORKTREE_PREFIX in validate_artifacts.py.
EXTENSION_ROOT = ROOT / "results" / "v2" / "extension"
CANONICAL_V2_RESULTS = ROOT / "results" / "v2"
HEAVY_ARTIFACT_ROOT = "results/v2/extension/"
CLAIM_PACK_ROOT = "docs/reports/extension-claim-pack/"
PLAN_ROOT = ROOT / "docs" / "reports" / "extension-plans"

# Fail-closed write confinement for tags / --out-dir (no .. / no escape).
WRITE_CONFINEMENT_ROOTS: tuple[Path, ...] = (
    ROOT / "docs" / "reports" / "extension-claim-pack",
    ROOT / "results" / "v2" / "extension",
)
_SAFE_TAG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
# Peak names accepted by build_atac_graph_v2.py: chrom:start-end
_PEAK_COLON_RE = re.compile(r"^[^:]+:\d+-\d+$")
# OpenProblems / Signac-style: chrom-start-end
_PEAK_HYPHEN_RE = re.compile(r"^(.+)-(\d+)-(\d+)$")

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
    "descartes_spleen": [
        DESKTOP_DATA
        / "datasets"
        / "extension_pilots"
        / "descartes_spleen"
        / "descartes_spleen_peaks.h5ad",
        DESKTOP_DATA / "datasets" / "extension_pilots" / "descartes_spleen",
    ],
}

DESCARTES_PILOT_DIR = (
    DESKTOP_DATA / "datasets" / "extension_pilots" / "descartes_spleen"
)


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


def assert_safe_tag(tag: str, *, label: str = "tag") -> str:
    """Reject path separators and ``..`` in construct / baseline tags."""
    if not tag or ".." in tag or "/" in tag or "\\" in tag:
        raise ValueError(f"{label} rejects path separators or '..': {tag!r}")
    if not _SAFE_TAG_RE.match(tag):
        raise ValueError(f"{label} must match {_SAFE_TAG_RE.pattern}: {tag!r}")
    return tag


def assert_confined_write_path(path: Path | str, *, label: str = "path") -> Path:
    """Resolve ``path`` and require it under allowlisted write roots.

    Rejects ``..`` components and absolute/relative escapes outside
    ``docs/reports/extension-claim-pack/`` and ``results/v2/extension/``.
    """
    raw = Path(path)
    if ".." in raw.parts:
        raise ValueError(f"{label} rejects '..' path components: {path}")
    resolved = raw.resolve() if raw.is_absolute() else (ROOT / raw).resolve()
    for root in WRITE_CONFINEMENT_ROOTS:
        root_r = root.resolve()
        try:
            if resolved == root_r or resolved.is_relative_to(root_r):
                return resolved
        except (ValueError, OSError):
            continue
    allowed = ", ".join(
        r.relative_to(ROOT).as_posix() + "/" for r in WRITE_CONFINEMENT_ROOTS
    )
    raise ValueError(f"{label} outside allowlisted roots ({allowed}): {path}")


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


def peerj_support_lock_enabled(env: Mapping[str, str] | None = None) -> bool:
    e = env if env is not None else os.environ
    return str(e.get("SCREG_PEERJ_SUPPORT_LOCK", "")).strip() == "1"


def under_extension_overlay(path: Path | str) -> bool:
    """True if ``path`` resolves under ``results/v2/extension/``."""
    resolved = Path(path).resolve()
    root = EXTENSION_ROOT.resolve()
    try:
        return resolved == root or resolved.is_relative_to(root)
    except (ValueError, OSError):
        return False


def is_canonical_v2_results_path(path: Path | str) -> bool:
    """True for ``results/v2/<file>`` but not ``results/v2/extension/...``."""
    resolved = Path(path).resolve()
    canon = CANONICAL_V2_RESULTS.resolve()
    try:
        if resolved == canon:
            return True
        if not resolved.is_relative_to(canon):
            return False
    except (ValueError, OSError):
        return False
    return not under_extension_overlay(resolved)


def normalize_peak_name(name: str) -> str:
    """Normalize peak ids to ``chrom:start-end`` for ``build_atac_graph_v2``."""
    s = str(name).strip()
    if _PEAK_COLON_RE.match(s):
        return s
    m = _PEAK_HYPHEN_RE.match(s)
    if m:
        return f"{m.group(1)}:{m.group(2)}-{m.group(3)}"
    return s


def peak_name_is_builder_ready(name: str) -> bool:
    return bool(_PEAK_COLON_RE.match(str(name).strip()))


def resolve_builder_out_dir(
    *,
    extension_out: str | None = None,
    peerj_lock: bool | None = None,
    env: Mapping[str, str] | None = None,
    create: bool = True,
) -> Path:
    """Resolve NPZ/meta output directory for ``build_atac_graph_v2``.

    - If ``SCREG_EXTENSION_OUT`` / ``extension_out`` is set → write there
      (mkdir parents).
    - If ``SCREG_PEERJ_SUPPORT_LOCK=1`` → refuse canonical ``results/v2/``
      writes; require an overlay path under ``results/v2/extension/``.
    - If both unset → default ``results/v2/`` (backward compatible).
    """
    e = env if env is not None else os.environ
    if extension_out is None:
        raw_ext = str(e.get("SCREG_EXTENSION_OUT", "")).strip()
        extension_out = raw_ext or None
    if peerj_lock is None:
        peerj_lock = peerj_support_lock_enabled(e)

    default = CANONICAL_V2_RESULTS

    if extension_out:
        raw = Path(extension_out)
        if ".." in raw.parts:
            raise ValueError(f"SCREG_EXTENSION_OUT rejects '..': {extension_out}")
        out = raw.resolve() if raw.is_absolute() else (ROOT / raw).resolve()
        if peerj_lock:
            if not under_extension_overlay(out):
                raise ValueError(
                    "SCREG_PEERJ_SUPPORT_LOCK=1 requires SCREG_EXTENSION_OUT under "
                    f"{HEAVY_ARTIFACT_ROOT} (refusing {redact_path(out) or out})"
                )
            # Reuse confinement helper for absolute/relative escape checks.
            out = assert_confined_write_path(out, label="SCREG_EXTENSION_OUT")
        if create:
            out.mkdir(parents=True, exist_ok=True)
        return out

    if peerj_lock:
        raise ValueError(
            "SCREG_PEERJ_SUPPORT_LOCK=1 requires SCREG_EXTENSION_OUT under "
            f"{HEAVY_ARTIFACT_ROOT} (refusing canonical results/v2/ PeerJ freeze paths)"
        )

    if create:
        default.mkdir(parents=True, exist_ok=True)
    return default.resolve()


def builder_env_command(
    *,
    tag: str,
    atac_file: str,
    extension_out: str,
    meta_file: str | None = None,
    extra_env: Mapping[str, str] | None = None,
) -> str:
    """Shell one-liner to build G_ATAC into the extension overlay (no fetch)."""
    tag = assert_safe_tag(tag, label="TAG")
    assert_confined_write_path(extension_out, label="SCREG_EXTENSION_OUT")
    parts = [
        f"TAG={tag}",
        f"ATAC_FILE={atac_file}",
        f"SCREG_EXTENSION_OUT={extension_out}",
        "SCREG_PEERJ_SUPPORT_LOCK=1",
    ]
    if meta_file:
        parts.append(f"META_FILE={meta_file}")
    if extra_env:
        for k, v in extra_env.items():
            parts.append(f"{k}={v}")
    parts.append("python src/v2/build_atac_graph_v2.py")
    return " ".join(parts)
