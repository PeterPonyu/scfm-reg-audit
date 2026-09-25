#!/usr/bin/env python
"""Shared path resolution for scReg-Eval scripts.

Precedence: environment variable -> project-local ./data -> error with the
variable name. No personal machine paths anywhere in the source tree.

Environment variables (see ENVIRONMENT.example at the project root):
  SCREG_DATA_ROOT   external data lake (models, ATAC datasets)
  SCREG_MODEL_ROOT  model weights root (default: $SCREG_DATA_ROOT/models)
  SCFM_BRAIN_ATAC   brain snATAC H5AD (default: $SCREG_DATA_ROOT/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad)
"""
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

_BRAIN_ATAC_REL = os.path.join(
    "datasets", "ATAC_data", "GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")


def project_root():
    return PROJECT_ROOT


def data_root():
    """External data root: SCREG_DATA_ROOT or <project>/data."""
    return Path(os.environ.get(
        "SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data")))


def model_root():
    return Path(os.environ.get("SCREG_MODEL_ROOT", data_root() / "models"))


def brain_atac_path():
    return os.environ.get("SCFM_BRAIN_ATAC", str(data_root() / _BRAIN_ATAC_REL))


def proxy_metadata_path(root=None):
    """Require project/configured metadata; only META_FILE=none selects pooled mode."""
    base = Path(root) if root is not None else PROJECT_ROOT
    configured = os.environ.get("META_FILE")
    if configured == "none":
        return None
    path = Path(configured) if configured is not None else base / "data/annotation/atac_cell_meta.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(
            f"ATAC cell-type metadata missing: {path}. Set META_FILE to the metadata "
            "CSV.gz, or explicitly set META_FILE=none for a pooled proxy; "
            "pooled and cell-type proxies are different analyses.")
    return str(path)
