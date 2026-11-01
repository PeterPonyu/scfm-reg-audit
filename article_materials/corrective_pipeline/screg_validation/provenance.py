"""Source candidates remain candidates until accession and transform evidence agree."""
from __future__ import annotations

import numpy as np
from scipy import sparse

from .common import ContractError, digest, identities, unique


def donor_join(cells, rows, *, cell_key="cell", donor_key="donor"):
    cells = unique(cells, "RNA cell key")
    lookup = {}
    for row in rows:
        cell = identities([row.get(cell_key)], "metadata cell key")[0]
        donor = identities([row.get(donor_key)], "donor")[0]
        if cell in lookup:
            raise ContractError("ambiguous barcode/donor join; composite sample+barcode keys required")
        lookup[cell] = donor
    missing = set(cells) - set(lookup)
    if missing:
        raise ContractError(f"unmatched cells in exact donor join: {len(missing)}")
    return [lookup[cell] for cell in cells]


def compare_candidate(left_cells, left_genes, left, right_cells, right_genes, right, *, atol=0., rtol=0.):
    lc, lg = unique(left_cells, "left cell"), unique(left_genes, "left gene")
    rc, rg = unique(right_cells, "right cell"), unique(right_genes, "right gene")
    left, right = sparse.csr_matrix(left), sparse.csr_matrix(right)
    if left.shape != (len(lc), len(lg)) or right.shape != (len(rc), len(rg)):
        raise ContractError("source fingerprint shape mismatch")
    ri, rj = {s: i for i, s in enumerate(rc)}, {s: i for i, s in enumerate(rg)}
    if not set(lc) <= set(rc) or not set(lg) <= set(rg):
        return {"status": "candidate_rejected", "all_cells_matched": set(lc) <= set(rc),
                "all_genes_matched": set(lg) <= set(rg), "provenance_closed": False}
    aligned = right[[ri[s] for s in lc]][:, [rj[s] for s in lg]]
    difference = (left - aligned).tocoo()
    if not np.isfinite(left.data).all() or not np.isfinite(aligned.data).all():
        raise ContractError("nonfinite source fingerprint")
    maxerr = float(np.max(np.abs(difference.data))) if difference.nnz else 0.
    scale = aligned[difference.row, difference.col].A1 if difference.nnz else np.array([])
    matched = bool(np.all(np.abs(difference.data) <= atol + rtol * np.abs(scale)))
    return {"status": "matrix_candidate_match" if matched else "candidate_rejected",
            "all_cells_matched": True, "all_genes_matched": True, "max_abs_error": maxerr,
            "n_cells": len(lc), "n_genes": len(lg), "atol": atol, "rtol": rtol,
            "matched_identity_sha256": digest([lc, lg]), "provenance_closed": False,
            "scope": "numeric match requires recorded accession, preprocessing chain and donor map before source closure"}


def validate_chain(chain):
    steps = chain.get("steps", [])
    if not chain.get("accession") or not chain.get("source_evidence_sha256") or not steps:
        raise ContractError("accession, source evidence and preprocessing steps required")
    for i, step in enumerate(steps):
        for field in ("input_sha256", "output_sha256", "command"):
            if not step.get(field):
                raise ContractError(f"missing provenance step {field}")
        for field in ("input_sha256", "output_sha256"):
            if len(step[field]) != 64 or any(c not in "0123456789abcdef" for c in step[field]):
                raise ContractError("invalid provenance hash")
        if i and steps[i-1]["output_sha256"] != step["input_sha256"]:
            raise ContractError("broken preprocessing hash chain")
    return {"chain_sha256": digest(chain), "chain_links_valid": True,
            "scope": "structural validation; source evidence and command replay require separate verification"}


def training_overlap(training, evaluation):
    """Missing listing is unknown; exclusion needs complete checkpoint-bound identities."""
    train = set(identities(training.get("dataset_ids", []), "training dataset"))
    evaluated = set(identities(evaluation.get("dataset_ids", []), "evaluation dataset"))
    if not evaluated:
        raise ContractError("evaluation dataset identities missing")
    overlap = sorted(train & evaluated)
    bound = (training.get("checkpoint_sha256") == evaluation.get("checkpoint_sha256")
             and bool(training.get("checkpoint_sha256")) and bool(training.get("source_evidence_sha256")))
    if overlap and bound:
        status = "documented_dataset_overlap"
    elif (bound and training.get("complete_dataset_inventory") is True
          and evaluation.get("identity_crosswalk_complete") is True):
        status = "excluded_at_documented_dataset_identity_resolution"
    else:
        status = "unknown"
    return {"status": status, "matching_dataset_ids": overlap, "checkpoint_bound": bool(bound),
            "never_pretrained_claim_allowed": False,
            "scope": "dataset-identity inventory only; donor/cell duplicates and undocumented aliases require their own audit"}
