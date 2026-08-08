#!/usr/bin/env python
"""Generate a provenance-bound pooled PBMC UCE gene graph."""
import hashlib
import json
import os
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata, spearmanr

import fm_readout as fr
import fm_readout_uce as fuce
import pbmc_cache


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data/manifest/shared_genes.v2.json"
RNA = ROOT / "data/multiome/pbmc10k_rna.h5ad"
G_ATAC = ROOT / "results/v2/G_ATAC_v2_PBMC10k.npz"
OUT = ROOT / "results/v2"
DATA_ROOT = Path(os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data")))
CHECKPOINT = DATA_ROOT / "models/UCE/4layer_model.torch"
ESM2_HUMAN = ROOT / "data/uce/human_esm2.pt"
SELECTION_SEED = 20260713
POOL_CAP = int(os.environ.get("UCE_POOL_CAP", "4000"))
BATCH = int(os.environ.get("UCE_BATCH", "4"))
CACHE_PATH = Path(os.environ.get("UCE_CACHE_PATH", OUT / "pbmc_uce_pooled_v2.npz"))
RESULT_PATH = Path(os.environ.get(
    "UCE_RESULT_PATH", str(CACHE_PATH.with_suffix(".json"))))
EXPECTED_CHECKPOINT_SHA256 = "acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560"
EXPECTED_ESM2_SHA256 = "a210e1cc7901513999b2bca3836ba9e2f203cd008be4e9a9d6412a2267de9748"
CO_NORMALIZATION_VERSION = "cp10k_log1p_v1"


def normalized_log_counts(matrix):
    """Return CP10k log1p values using the shared RNA normalization contract."""
    x = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    totals = np.asarray(x.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    cp10k = x.multiply(1e4 / totals[:, None]).tocsr()
    cp10k.data = np.log1p(cp10k.data)
    return cp10k


def log(*values):
    print(f"[{time.strftime('%H:%M:%S')}]", *values, flush=True)


def _atomic_json(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name, suffix=".json", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(document, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_legacy_uce_cache(path, cell_ids, genes, manifest_sha, pool_cap, rna_sha256,
                          checkpoint_sha256, esm2_sha256):
    """Load pre-contract cache only to migrate its UCE embedding, never its co graph."""
    with np.load(path, allow_pickle=False) as cache:
        metadata_matches = (
            np.array_equal(cache["cell_ids"], cell_ids)
            and [str(g) for g in cache["genes"]] == list(genes)
            and str(cache["manifest_sha"].item()) == manifest_sha
            and int(cache["selection_seed"].item()) == SELECTION_SEED
            and int(cache["pool_cap"].item()) == pool_cap
            and str(cache["rna_sha256"].item()) == rna_sha256
            and str(cache["checkpoint_sha256"].item()) == checkpoint_sha256
            and str(cache["esm2_sha256"].item()) == esm2_sha256
        )
        if not metadata_matches:
            raise ValueError("PBMC UCE legacy cache provenance mismatch")
        uce = cache["uce"].copy()
        covered = int(cache["covered"].item())
    shape = (len(genes), len(genes))
    if uce.shape != shape or not np.isfinite(uce).all():
        raise ValueError("PBMC UCE legacy cache graph shape mismatch")
    return uce, covered


def load_uce_cache(path, cell_ids, genes, manifest_sha, pool_cap, rna_sha256,
                   checkpoint_sha256, esm2_sha256):
    if not Path(path).exists():
        return None
    with np.load(path, allow_pickle=False) as cache:
        metadata_matches = (
            np.array_equal(cache["cell_ids"], cell_ids)
            and [str(g) for g in cache["genes"]] == list(genes)
            and str(cache["manifest_sha"].item()) == manifest_sha
            and int(cache["selection_seed"].item()) == SELECTION_SEED
            and int(cache["pool_cap"].item()) == pool_cap
            and str(cache["rna_sha256"].item()) == rna_sha256
            and str(cache["checkpoint_sha256"].item()) == checkpoint_sha256
            and str(cache["esm2_sha256"].item()) == esm2_sha256
            and str(cache["co_normalization_version"].item()) == CO_NORMALIZATION_VERSION
        )
        if not metadata_matches:
            raise ValueError("PBMC UCE cache provenance mismatch")
        co = cache["co"].copy()
        uce = cache["uce"].copy()
        covered = int(cache["covered"].item())
    shape = (len(genes), len(genes))
    if co.shape != shape or uce.shape != shape:
        raise ValueError("PBMC UCE cache graph shape mismatch")
    if not np.isfinite(co).all() or not np.isfinite(uce).all():
        raise ValueError("PBMC UCE cache contains non-finite values")
    return co, uce, covered


def write_uce_cache(path, co, uce, covered, cell_ids, genes, manifest_sha,
                    pool_cap, rna_sha256, checkpoint_sha256, esm2_sha256):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name, suffix=".npz", dir=path.parent)
    os.close(fd)
    try:
        np.savez(
            temporary,
            co=co,
            uce=uce,
            covered=np.asarray(covered),
            cell_ids=np.asarray(cell_ids),
            genes=np.asarray(genes),
            manifest_sha=np.asarray(manifest_sha),
            selection_seed=np.asarray(SELECTION_SEED),
            pool_cap=np.asarray(pool_cap),
            rna_sha256=np.asarray(rna_sha256),
            checkpoint_sha256=np.asarray(checkpoint_sha256),
            esm2_sha256=np.asarray(esm2_sha256),
            co_normalization_version=np.asarray(CO_NORMALIZATION_VERSION),
        )
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def partial_spearman(x, y, control):
    rx, ry, rc = rankdata(x), rankdata(y), rankdata(control)
    design = np.column_stack([np.ones(len(rc)), rc])
    bx = np.linalg.lstsq(design, rx, rcond=None)[0]
    by = np.linalg.lstsq(design, ry, rcond=None)[0]
    return float(np.corrcoef(rx - design @ bx, ry - design @ by)[0, 1])


def main():
    manifest = json.loads(MANIFEST.read_text())
    genes = manifest["genes"]
    manifest_sha = manifest["sha256"]
    observed_manifest_sha = hashlib.sha256("\n".join(genes).encode()).hexdigest()
    if observed_manifest_sha != manifest_sha:
        raise ValueError("manifest self-hash mismatch")

    checkpoint_sha = pbmc_cache.sha256_file(CHECKPOINT)
    esm2_sha = pbmc_cache.sha256_file(ESM2_HUMAN)
    if checkpoint_sha != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError("UCE checkpoint SHA-256 mismatch")
    if esm2_sha != EXPECTED_ESM2_SHA256:
        raise ValueError("UCE human ESM2 SHA-256 mismatch")

    rna_sha = pbmc_cache.sha256_file(RNA)
    adata = ad.read_h5ad(RNA)
    cell_ids = pbmc_cache.select_pool_cell_ids(adata.n_obs, POOL_CAP, SELECTION_SEED)
    try:
        cached = load_uce_cache(
            CACHE_PATH, cell_ids, genes, manifest_sha, POOL_CAP, rna_sha,
            checkpoint_sha, esm2_sha,
        )
        legacy = None
    except KeyError:
        cached = None
        legacy = load_legacy_uce_cache(
            CACHE_PATH, cell_ids, genes, manifest_sha, POOL_CAP, rna_sha,
            checkpoint_sha, esm2_sha,
        )
        log("legacy PBMC UCE cache lacks normalization metadata; migrating UCE embedding only")

    if cached is None:
        matrix = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
        selected = matrix[cell_ids]
        if selected.data.size and not np.equal(selected.data, np.floor(selected.data)).all():
            raise ValueError("PBMC UCE input must contain raw integer counts")

        if legacy is None:
            readout = fuce.UCEReadout(genes, batch=BATCH, seed=SELECTION_SEED)
            symbols = {str(symbol): index for index, symbol in enumerate(adata.var_names)}
            universe_columns = np.array([symbols.get(gene, -1) for gene in readout.uni_syms])
            present = universe_columns >= 0
            xraw = matrix[:, np.where(present, universe_columns, 0)].tocsr()
            xraw = xraw.multiply(present[None, :]).tocsr()
            log(f"UCE over {len(cell_ids)} cells; universe coverage {int(present.sum())}/{len(present)}")
            embeddings = readout.gene_embed(xraw, cell_ids, universe_columns)
            covered = int(np.sum(np.abs(embeddings).sum(axis=1) > 0))
            uce = fr.FMReadout.cos_graph(embeddings)
        else:
            symbols = {str(symbol): index for index, symbol in enumerate(adata.var_names)}
            uce, covered = legacy

        manifest_columns = np.array([symbols.get(gene, -1) for gene in genes])
        manifest_present = manifest_columns >= 0
        xlog = normalized_log_counts(
            matrix[cell_ids][:, np.where(manifest_present, manifest_columns, 0)])
        co = np.zeros((len(genes), len(genes)), dtype=np.float32)
        present_graph = fr.gene_coexp(xlog[:, manifest_present].toarray())
        co[np.ix_(manifest_present, manifest_present)] = present_graph
        write_uce_cache(
            CACHE_PATH, co, uce, covered, cell_ids, genes, manifest_sha,
            POOL_CAP, rna_sha, checkpoint_sha, esm2_sha,
        )
        log(f"saved cache {CACHE_PATH}")
    else:
        co, uce, covered = cached
        log(f"validated and reused cache {CACHE_PATH}")

    with np.load(G_ATAC, allow_pickle=False) as atac_cache:
        types = [str(value) for value in atac_cache["types"]]
        tf_rows = atac_cache["tf_rows"].copy()
        atac = np.mean([atac_cache[f"G_{cell_type}"] for cell_type in types], axis=0)
    ii = np.repeat(tf_rows, len(genes))
    jj = np.tile(np.arange(len(genes)), len(tf_rows))
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    atac_values, co_values, uce_values = atac[ii, jj], co[ii, jj], uce[ii, jj]
    metrics = {
        "coexp_vs_atac": float(spearmanr(co_values, atac_values).statistic),
        "uce_vs_atac": float(spearmanr(uce_values, atac_values).statistic),
        "uce_vs_coexp": float(spearmanr(uce_values, co_values).statistic),
        "uce_partial_given_coexp": partial_spearman(uce_values, atac_values, co_values),
    }
    result = {
        "schema_version": 2,
        "fm": "UCE_4layer",
        "tissue": "PBMC10k_multiome_paired",
        "n_cells": int(len(cell_ids)),
        "n_cells_coexp": int(len(cell_ids)),
        "pool_cap": POOL_CAP,
        "selection_seed": SELECTION_SEED,
        "manifest_genes_covered": covered,
        "n_pairs": int(len(ii)),
        "co_normalization_version": CO_NORMALIZATION_VERSION,
        "metrics": {key: round(value, 6) for key, value in metrics.items()},
        "cache": str(CACHE_PATH),
        "provenance": {
            "manifest_sha256": manifest_sha,
            "rna_sha256": rna_sha,
            "checkpoint_sha256": checkpoint_sha,
            "esm2_sha256": esm2_sha,
            "cache_sha256": pbmc_cache.sha256_file(CACHE_PATH),
        },
    }
    _atomic_json(RESULT_PATH, result)
    log(json.dumps(result["metrics"], sort_keys=True))
    log(f"saved result {RESULT_PATH}")


if __name__ == "__main__":
    main()
