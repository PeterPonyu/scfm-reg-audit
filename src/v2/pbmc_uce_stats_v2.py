#!/usr/bin/env python
"""Append the matched-control PBMC UCE row to accepted pooled statistics."""
import copy
import json
import os
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np

import fixed_panel_audit as fpa
import pbmc_cache
import run_fixed_panel_audit as drv


N_PERM_MANTEL = int(os.environ.get("N_PERM_POOLED_MANTEL", "999"))
N_PERM_DEG = int(os.environ.get("N_PERM_POOLED_DEG", "999"))
BASE_PATH = Path(os.environ.get(
    "PBMC_UCE_STATS_BASE", f"{fpa.OUT}/pbmc_scgpt_stats_v2.json"))
UCE_PATH = Path(os.environ.get(
    "PBMC_UCE_CACHE", f"{fpa.OUT}/pbmc_uce_pooled_v2.npz"))
OUT_PATH = Path(os.environ.get(
    "PBMC_UCE_STATS_OUT", f"{fpa.OUT}/pbmc_uce_stats_v2.json"))
SELECTION_SEED = 20260713
POOL_CAP = 4000


def log(*values):
    print(f"[{time.strftime('%H:%M:%S')}]", *values, flush=True)


def write_json_atomic(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name, suffix=".json", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(document, handle, indent=2, allow_nan=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_uce_graph(path):
    genes, _detection, manifest_sha = fpa.load_manifest()
    rna_path = Path(fpa.ROOT) / "data/multiome/pbmc10k_rna.h5ad"
    rna = ad.read_h5ad(rna_path, backed="r")
    expected_ids = pbmc_cache.select_pool_cell_ids(rna.n_obs, POOL_CAP, SELECTION_SEED)
    rna.file.close()
    expected_rna_sha = pbmc_cache.sha256_file(rna_path)
    with np.load(path, allow_pickle=False) as cache:
        if [str(gene) for gene in cache["genes"]] != genes:
            raise ValueError("PBMC UCE graph gene order mismatch")
        if str(cache["manifest_sha"].item()) != manifest_sha:
            raise ValueError("PBMC UCE graph manifest mismatch")
        if int(cache["selection_seed"].item()) != SELECTION_SEED:
            raise ValueError("PBMC UCE graph selection seed mismatch")
        if int(cache["pool_cap"].item()) != POOL_CAP:
            raise ValueError("PBMC UCE graph pool cap mismatch")
        if not np.array_equal(cache["cell_ids"], expected_ids):
            raise ValueError("PBMC UCE graph selected cells mismatch")
        if str(cache["rna_sha256"].item()) != expected_rna_sha:
            raise ValueError("PBMC UCE graph RNA input mismatch")
        co, uce = cache["co"].copy(), cache["uce"].copy()
    expected_shape = (len(genes), len(genes))
    if co.shape != expected_shape or uce.shape != expected_shape:
        raise ValueError("PBMC UCE graph shape mismatch")
    if not np.isfinite(co).all() or not np.isfinite(uce).all():
        raise ValueError("PBMC UCE graph contains non-finite values")
    return co, uce


def main():
    base = json.loads(BASE_PATH.read_text())
    if base.get("n_perm_mantel") != N_PERM_MANTEL:
        raise ValueError("base Mantel permutation count mismatch")
    if base.get("n_perm_degree") != N_PERM_DEG:
        raise ValueError("base degree permutation count mismatch")
    result = copy.deepcopy(base["pooled_pbmc"])
    existing = {row.get("model_label") for row in result["primary_family"]["rows"]}
    if "UCE_encoder" in existing:
        raise ValueError("base statistics already contain UCE_encoder")

    g_atac, _co, _models, tf_rows, types = drv.load_pooled_pbmc()
    uce_co, uce_graph = load_uce_graph(UCE_PATH)
    atac_path = f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"
    drv.append_independent_control_model(
        result, atac_path, "pbmc", g_atac, uce_co, uce_graph, tf_rows,
        "UCE_encoder", N_PERM_MANTEL, N_PERM_DEG,
        np.random.SeedSequence([drv.SEED_ROOT, 20260729, 2]), str(UCE_PATH),
        model_key="uce", control_key="co",
    )
    document = {
        "schema_version": 1,
        "design": "PBMC pooled fixed-panel update with matched-control scGPT and UCE",
        "seed_root": drv.SEED_ROOT,
        "n_perm_mantel": N_PERM_MANTEL,
        "n_perm_degree": N_PERM_DEG,
        "base_statistics": str(BASE_PATH),
        "base_statistics_sha256": fpa.sha256_file(BASE_PATH),
        "types": types,
        "pooled_pbmc": result,
    }
    write_json_atomic(OUT_PATH, document)
    log(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
