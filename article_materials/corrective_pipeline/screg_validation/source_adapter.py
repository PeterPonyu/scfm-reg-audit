"""Bounded candidate-source comparison using explicit H5AD matrices and scales."""
from __future__ import annotations

import hashlib

import h5py
import numpy as np

from .common import ContractError, sha256, unique
from .h5ad import counts_normalization, index, read_rows
from .provenance import compare_candidate


def compare_h5ad_sources(left, right, *, max_cells, max_genes, seed, atol, rtol):
    """Deterministic fingerprint, not a full matrix/source provenance certification."""
    if max_cells < 1 or max_genes < 1:
        raise ContractError("source fingerprint budget required")
    metadata = []
    for resource in (left, right):
        if resource.get("sha256") != sha256(resource["path"]):
            raise ContractError("source candidate hash drift")
        with h5py.File(resource["path"], "r") as f:
            metadata.append((index(f, "obs"), index(f, "var")))
    (lc, lg), (rc, rg) = metadata
    if not set(lc) <= set(rc) or not set(lg) <= set(rg):
        return {"status": "candidate_rejected", "full_cell_identity_subset": set(lc) <= set(rc),
                "full_gene_identity_subset": set(lg) <= set(rg), "provenance_closed": False,
                "reason": "barcode or feature universe differs; require an evidence-backed identity crosswalk, never suffix stripping"}
    chosen_cells = sorted(lc, key=lambda s: hashlib.sha256(f"{seed}:{s}".encode()).hexdigest())[:max_cells]
    chosen_genes = sorted(lg, key=lambda s: hashlib.sha256(f"{seed}:{s}".encode()).hexdigest())[:max_genes]
    matrices = []
    for resource, (cells, genes) in zip((left, right), metadata):
        ci, gi = {c: i for i, c in enumerate(cells)}, {g: i for i, g in enumerate(genes)}
        values = read_rows(resource["path"], resource["matrix"], [ci[c] for c in chosen_cells])
        _, log = counts_normalization(values, resource["scale"])
        matrices.append(log[:, [gi[g] for g in chosen_genes]])
    result = compare_candidate(chosen_cells, chosen_genes, matrices[0], chosen_cells, chosen_genes, matrices[1], atol=atol, rtol=rtol)
    result.update(full_cell_identity_subset=True, full_gene_identity_subset=True,
                  left_input_sha256=left["sha256"], right_input_sha256=right["sha256"],
                  comparison_scale="explicit log1p CP10k, normalized before selecting fingerprint genes",
                  scope="bounded matrix fingerprint plus complete identity subset; matching does not prove accession, donor identity or preprocessing provenance")
    return result
