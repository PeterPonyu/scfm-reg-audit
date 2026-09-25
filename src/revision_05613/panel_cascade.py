#!/usr/bin/env python
"""R2.1: filter cascade behind the frozen 1,200-gene / 446-TF panel.

Re-applies the rules of src/v2/freeze_gene_manifest.py read-only (nothing is
written to data/manifest) and counts genes and TFs after each step, then checks
that the capped panel reproduces the frozen manifest and its SHA-256.
"""
import json
import pickle
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import log  # noqa: E402

sys.path.insert(0, str(common.SRC / "v2"))
import motif_utils as mu  # noqa: E402

ROOT = Path(common.fpa.ROOT)
DATA_ROOT = Path(common.rfa.DATA_ROOT)
COORDS = ROOT / "data/annotation/gene_coords_hg38.tsv"
MEME = ROOT / "data/motifs/JASPAR2024_CORE_vertebrates.meme"
ATAC = DATA_ROOT / "datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
RNA = DATA_ROOT / "datasets/extra_preprocessed/ad_hm_prepped.h5ad"
GF = DATA_ROOT / "models/Geneformer/geneformer"
CK = DATA_ROOT / "models/scGPT-human"
PROMOTER, DETECT_FRAC, N_MAX = 2000, 0.01, 1200


def main():
    man = json.loads(Path(common.fpa.MANI).read_text())
    gtok = pickle.load(open(GF / "token_dictionary_gc104M.pkl", "rb"))
    gmed = pickle.load(open(GF / "gene_median_dictionary_gc104M.pkl", "rb"))
    gn2i = pickle.load(open(GF / "gene_name_id_dict_gc104M.pkl", "rb"))
    svoc = json.load(open(CK / "vocab.json"))

    def gf_ok(s):
        e = gn2i.get(s)
        return e is not None and e in gtok and e in gmed

    mot = mu.parse_meme(str(MEME))
    TF = {t for d in mot.values() for t in mu.tf_symbols(d["name"])}

    genes, n_records = {}, 0
    for ln in open(COORDS):
        chrom, s, e, strand, name = ln.rstrip("\n").split("\t")
        n_records += 1
        if name in genes:
            continue
        s, e = int(s), int(e)
        lo = s - PROMOTER if strand == "+" else s
        hi = e if strand == "+" else e + PROMOTER
        genes[name] = (chrom, lo, hi)

    Av = ad.read_h5ad(ATAC, backed="r")
    names = [str(p) for p in Av.var_names]
    n_peaks = len(names)
    pchr = np.array([p.split(":")[0] for p in names])
    pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in names])
    pmid = (pse[:, 0] + pse[:, 1]) // 2
    by_chr = {}
    for i, c in enumerate(pchr):
        by_chr.setdefault(c, []).append(i)
    by_chr = {c: np.array(v) for c, v in by_chr.items()}
    has_peak = set()
    for name, (chrom, lo, hi) in genes.items():
        pis = by_chr.get(chrom)
        if pis is not None and np.any((pmid[pis] >= lo) & (pmid[pis] <= hi)):
            has_peak.add(name)

    R = ad.read_h5ad(RNA)
    X = R.X.tocsc() if sp.issparse(R.X) else sp.csr_matrix(R.X).tocsc()
    ncell = X.shape[0]
    det = {str(s): float(X[:, j].getnnz()) / ncell for j, s in enumerate(R.var_names)}

    steps = []

    def record(label, pool):
        pool = list(pool)
        steps.append({"step": label, "genes": len(pool), "tfs": sum(g in TF for g in pool)})

    pool = list(genes)
    record("GENCODE v44 basic symbols with coordinates (first record per symbol)", pool)
    pool = [g for g in pool if gf_ok(g)]
    record("Geneformer token with a median (gc104M dictionaries)", pool)
    pool = [g for g in pool if g in svoc]
    record("in the scGPT whole-human vocabulary", pool)
    pool = [g for g in pool if g in has_peak]
    record(">=1 GSE174367 snATAC peak midpoint in [TSS-2 kb, TES] (strand-aware)", pool)
    pool = [g for g in pool if det.get(g, 0.0) >= DETECT_FRAC]
    record(f"RNA detection >= {DETECT_FRAC:.0%} of cells in the brain RNA reference", pool)
    universe = pool

    uni_tf = sorted(g for g in universe if g in TF)
    uni_non = sorted((g for g in universe if g not in TF), key=lambda g: -det[g])
    chosen = sorted(set(uni_tf + uni_non[: max(0, N_MAX - len(uni_tf))]))
    record(f"cap to {N_MAX}: all TFs kept, remaining slots by RNA detection (descending)", chosen)

    import hashlib
    sha = hashlib.sha256("\n".join(chosen).encode()).hexdigest()
    reproduced = chosen == man["genes"] and sha == man["sha256"]
    non_tf_det = sorted((det[g] for g in chosen if g not in TF))
    result = {
        "schema_version": 1,
        "analysis": "panel_filter_cascade",
        "source_rules": "src/v2/freeze_gene_manifest.py (re-applied read-only)",
        "coordinate_records": n_records,
        "jaspar_tf_symbols": len(TF),
        "brain_atac_peaks": n_peaks,
        "rna_cells": int(ncell),
        "steps": steps,
        "universe_genes": len(universe),
        "universe_tfs": len(uni_tf),
        "non_tf_slots": N_MAX - len(uni_tf),
        "non_tf_detection_min_filled": float(non_tf_det[0]),
        "manifest_sha256": man["sha256"],
        "reproduced_sha256": sha,
        "reproduces_frozen_manifest": bool(reproduced),
        "note": ("The coordinate file covers all GENCODE v44 basic biotypes, not only "
                 "protein-coding genes; the FM vocabularies do most of the narrowing."),
    }
    common.write_json(common.OUT_DIR / "panel_cascade.json", result)
    for s in steps:
        log(f"{s['genes']:>7} genes {s['tfs']:>5} TFs  {s['step']}")
    log(f"reproduces frozen manifest: {reproduced} (sha {sha[:8]})")
    if not reproduced:
        sys.exit(2)


if __name__ == "__main__":
    main()
