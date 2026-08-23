import hashlib
import importlib.metadata
import json
import os
import platform
import tempfile

import numpy as np


def write_preflight_report(path, metrics, reference_path, manifest_sha):
    import anndata
    import pandas
    import scipy

    with open(reference_path, "rb") as fh:
        reference_sha = hashlib.sha256(fh.read()).hexdigest()
    report = {
        "status": "passed",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "anndata": importlib.metadata.version("anndata"),
        "manifest_sha": manifest_sha,
        "reference_metrics_sha256": reference_sha,
        "metrics": metrics,
    }
    parent = os.path.dirname(path) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=os.path.basename(path), suffix=".json", dir=parent)
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
            fh.write("\n")
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def verify_reference_metrics(current, path):
    with open(path) as fh:
        reference = json.load(fh)
    keys = ["n_pairs", "n_tf", "tissue", "types"]
    keys.extend(sorted(key for key in current if key.startswith(("embed__", "attn__"))))
    mismatches = {
        key: {"expected": reference.get(key), "observed": current.get(key)}
        for key in keys if reference.get(key) != current.get(key)
    }
    if mismatches:
        raise RuntimeError(f"PBMC reference metric drift: {mismatches}")


def load_confound_cache(path, genes, manifest_sha):
    with np.load(path, allow_pickle=False) as cache:
        cached_genes = [str(g) for g in cache["genes"]]
        cached_sha = str(cache["manifest_sha"].item())
        if cached_genes != list(genes) or cached_sha != manifest_sha:
            raise ValueError("PBMC confound cache manifest mismatch")
        vectors = tuple(cache[key].copy() for key in ("peakcount", "genelen", "detv", "gc"))
    if any(vector.shape != (len(genes),) for vector in vectors):
        raise ValueError("PBMC confound cache vector shape mismatch")
    if any(not np.isfinite(vector).all() for vector in vectors):
        raise ValueError("PBMC confound cache contains non-finite values")
    return vectors


def select_pool_cell_ids(n_cells, pool_cap, seed):
    cell_ids = np.arange(n_cells)
    if n_cells > pool_cap:
        cell_ids = np.random.default_rng(seed).choice(
            cell_ids, size=pool_cap, replace=False)
    return cell_ids


def sha256_file(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_scgpt_cache(path, cell_ids, genes, manifest_sha, selection_seed, pool_cap, rna_sha256):
    if not os.path.exists(path):
        return None
    with np.load(path, allow_pickle=False) as cache:
        cached_ids = cache["cell_ids"]
        cached_genes = [str(g) for g in cache["genes"]]
        cached_sha = str(cache["manifest_sha"].item())
        cached_seed = int(cache["selection_seed"].item())
        cached_cap = int(cache["pool_cap"].item())
        cached_rna_sha = str(cache["rna_sha256"].item())
        if (not np.array_equal(cached_ids, cell_ids) or cached_genes != list(genes)
                or cached_sha != manifest_sha or cached_seed != selection_seed
                or cached_cap != pool_cap or cached_rna_sha != rna_sha256):
            raise ValueError("PBMC scGPT cache provenance mismatch")
        co, sg = cache["co"].copy(), cache["sg"].copy()
    expected_shape = (len(genes), len(genes))
    if co.shape != expected_shape or sg.shape != expected_shape:
        raise ValueError("PBMC scGPT cache graph shape mismatch")
    if not np.isfinite(co).all() or not np.isfinite(sg).all():
        raise ValueError("PBMC scGPT cache contains non-finite values")
    return co, sg


def write_scgpt_cache(path, co, sg, cell_ids, genes, manifest_sha, selection_seed, pool_cap,
                       rna_sha256):
    parent = os.path.dirname(path) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=os.path.basename(path), suffix=".npz", dir=parent)
    os.close(fd)
    try:
        np.savez(
            tmp_path,
            co=co,
            sg=sg,
            cell_ids=np.asarray(cell_ids),
            genes=np.asarray(genes),
            manifest_sha=np.asarray(manifest_sha),
            selection_seed=np.asarray(selection_seed),
            pool_cap=np.asarray(pool_cap),
            rna_sha256=np.asarray(rna_sha256),
        )
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
