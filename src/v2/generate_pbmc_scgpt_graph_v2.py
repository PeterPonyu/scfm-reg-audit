#!/usr/bin/env python
import hashlib
import json
import os
import time

import numpy as np

import fm_readout as fr
import pbmc_cache


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANIFEST = f"{ROOT}/data/manifest/shared_genes.v2.json"
RNA = f"{ROOT}/data/multiome/pbmc10k_rna.h5ad"
OUT = f"{ROOT}/results/v2/pbmc_scgpt_pooled_v2.npz"
POOL_CAP = 4000
SELECTION_SEED = 20260713
BATCH = int(os.environ.get("FM_BATCH", "64"))


def log(*args):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, flush=True)


def main():
    manifest = json.load(open(MANIFEST))
    genes = manifest["genes"]
    manifest_sha = hashlib.sha256(("\n".join(genes)).encode()).hexdigest()
    if manifest_sha != manifest["sha256"]:
        raise ValueError("manifest hash mismatch")
    A, _, Xlog, _ = fr.load_norm(RNA)
    cell_ids = pbmc_cache.select_pool_cell_ids(A.shape[0], POOL_CAP, SELECTION_SEED)
    rna_sha = pbmc_cache.sha256_file(RNA)
    cached = pbmc_cache.load_scgpt_cache(
        OUT, cell_ids, genes, manifest_sha, SELECTION_SEED, POOL_CAP, rna_sha)
    if cached is not None:
        log(f"cache already valid: {OUT}")
        return
    gene_index = {str(symbol): index for index, symbol in enumerate(A.var_names)}
    present = np.array([index for index, gene in enumerate(genes) if gene in gene_index])
    source_columns = np.array([gene_index[genes[index]] for index in present])
    X = Xlog[cell_ids][:, source_columns].tocsr()
    log(f"scGPT graph: cells={len(cell_ids)} genes={len(present)} batch={BATCH}")
    co_small = fr.gene_coexp(X.toarray())
    readout = fr.FMReadout([genes[index] for index in present], batch=BATCH)
    embeddings = readout.scgpt(X, np.arange(X.shape[0]))
    sg_small = fr.FMReadout.cos_graph(embeddings)
    co = np.zeros((len(genes), len(genes)), dtype=np.float32)
    sg = np.zeros((len(genes), len(genes)), dtype=np.float32)
    co[np.ix_(present, present)] = co_small
    sg[np.ix_(present, present)] = sg_small
    pbmc_cache.write_scgpt_cache(
        OUT, co, sg, cell_ids, genes, manifest_sha, SELECTION_SEED,
        POOL_CAP, rna_sha)
    log(f"saved {OUT}")


if __name__ == "__main__":
    main()
