"""Freeze selection using discovery donors only; separate degree and evaluation."""
from __future__ import annotations

import hashlib
from itertools import combinations

import numpy as np
from scipy import sparse

from .common import ContractError, digest, identities, unique


def donor_split(donors, *, seed, discovery_count=3, reference_count=2, min_evaluation=3):
    donors = sorted(set(identities(donors, "donor")),
                    key=lambda s: hashlib.sha256(f"{seed}:{s}".encode()).hexdigest())
    if min(discovery_count, reference_count, min_evaluation) < 1 or len(donors) < discovery_count + reference_count + min_evaluation:
        raise ContractError("insufficient distinct donors for discovery/reference/evaluation")
    return {"discovery": donors[:discovery_count],
            "degree_reference": donors[discovery_count:discovery_count + reference_count],
            "evaluation": donors[discovery_count + reference_count:]}


def validate_split(split, observed_donors):
    if set(split) != {"discovery", "degree_reference", "evaluation"}:
        raise ContractError("three explicit donor partitions required")
    for name, donors in split.items():
        if not unique(donors, name + " donor"):
            raise ContractError("empty donor partition")
    for a, b in combinations(split, 2):
        if set(split[a]) & set(split[b]):
            raise ContractError("donor overlap across partitions")
    if set(sum(split.values(), [])) != set(identities(observed_donors, "observed donor")):
        raise ContractError("donor partition incomplete")


def select_panel(discovery_counts, genes, donor_ids, split, *, eligible_genes, tf_symbols,
                 original_genes, max_genes, min_genes, min_tfs, detection_floor, gene_disjoint):
    genes = unique(genes, "gene symbol")
    donors = identities(donor_ids, "discovery donor")
    x = sparse.csr_matrix(discovery_counts)
    if x.shape != (len(donors), len(genes)) or set(donors) != set(split["discovery"]):
        raise ContractError("panel selector accepts all and only discovery donors")
    if any(set(donors) & set(split[k]) for k in ("degree_reference", "evaluation")):
        raise ContractError("test/reference donor leaked into panel selection")
    if not 0 <= detection_floor <= 1 or not (0 < min_tfs <= min_genes <= max_genes):
        raise ContractError("invalid frozen selection thresholds")
    detected = (x > 0).astype(float)
    # Equal-donor detection avoids a donor with more cells choosing the panel.
    rate = np.mean([np.asarray(detected[np.array(donors) == d].mean(0)).ravel()
                    for d in split["discovery"]], axis=0)
    allowed = set(eligible_genes)
    if gene_disjoint:
        allowed -= set(original_genes)
    rates = dict(zip(genes, map(float, rate)))
    passed = [g for g in genes if g in allowed and rates[g] >= detection_floor]
    tf = set(tf_symbols)
    order = lambda g: (-rates[g], g)
    chosen = sorted([g for g in passed if g in tf], key=order)[:max_genes]
    chosen += sorted([g for g in passed if g not in tf], key=order)[:max_genes - len(chosen)]
    chosen = sorted(chosen)
    ntf = len(set(chosen) & tf)
    if len(chosen) < min_genes or ntf < min_tfs:
        raise ContractError("new panel fails preregistered coverage; do not relax thresholds after outcomes")
    body = {"genes": chosen, "tf_genes": sorted(set(chosen) & tf), "split": split,
            "detection": {g: rates[g] for g in chosen}, "gene_disjoint": gene_disjoint,
            "overlap_original": len(set(chosen) & set(original_genes)),
            "selection": "equal-donor RNA detection; TF priority; lexical ties; no evaluation RNA or outcome labels",
            "selection_independence": "new discovery rule and held-out donors"}
    return {**body, "panel_sha256": digest(body)}


def verify_panel(panel):
    body = {k: v for k, v in panel.items() if k != "panel_sha256"}
    if panel.get("panel_sha256") != digest(body):
        raise ContractError("frozen panel modified")
    unique(panel["genes"], "panel gene")
    validate_split(panel["split"], sum(panel["split"].values(), []))
    if not set(panel["tf_genes"]) <= set(panel["genes"]):
        raise ContractError("TFs outside frozen panel")
