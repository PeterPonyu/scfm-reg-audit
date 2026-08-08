import importlib.metadata
import json
import os
import platform

import numpy as np

import audit_utils as au

# Chunked file hashing is shared with every other cache writer.
sha256_file = au.sha256_file


def write_preflight_report(path, metrics, reference_path, manifest_sha):
    import anndata
    import pandas
    import scipy

    reference_sha = sha256_file(reference_path)
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
    au.write_json_atomic(path, report, allow_nan=True)


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


# ----------------------------- generic graph caches --------------------------
# Every model readout caches its (n_genes, n_genes) graphs next to the exact
# provenance that produced them: cell pool, gene panel, manifest hash, selection
# seed, pool cap, and input/checkpoint hashes. A cache whose provenance differs
# is never silently reused.
def check_provenance(cache, expectations, label):
    """Raise if any expected provenance field differs from the cached one.

    Arrays (cell ids) are compared element-wise, gene panels as string lists,
    and every other field by the expected value's own type.
    """
    for key, expected in expectations.items():
        cached = cache[key]
        if isinstance(expected, str):
            matches = str(cached.item()) == expected
        elif isinstance(expected, (int, np.integer)):
            matches = int(cached.item()) == int(expected)
        elif isinstance(expected, (list, tuple)):
            matches = [str(v) for v in cached] == [str(v) for v in expected]
        else:
            matches = np.array_equal(np.asarray(cached), np.asarray(expected))
        if not matches:
            raise ValueError(f"{label} cache provenance mismatch: {key}")


def load_graph_cache(path, graph_keys, expectations, label, n_genes, scalar_keys=()):
    """Load cached square graphs after a full provenance/shape/finiteness check.

    Returns None when the cache is absent, otherwise a dict of graphs plus any
    requested integer scalars (e.g. gene coverage).
    """
    if not os.path.exists(path):
        return None
    with np.load(path, allow_pickle=False) as cache:
        check_provenance(cache, expectations, label)
        out = {key: cache[key].copy() for key in graph_keys}
        out.update({key: int(cache[key].item()) for key in scalar_keys})
    expected_shape = (n_genes, n_genes)
    if any(out[key].shape != expected_shape for key in graph_keys):
        raise ValueError(f"{label} cache graph shape mismatch")
    if any(not np.isfinite(out[key]).all() for key in graph_keys):
        raise ValueError(f"{label} cache contains non-finite values")
    return out


def write_graph_cache(path, **arrays):
    """Atomically write graphs plus their provenance fields to an .npz."""
    au.write_npz_atomic(path, **{k: np.asarray(v) for k, v in arrays.items()})


def load_scgpt_cache(path, cell_ids, genes, manifest_sha, selection_seed, pool_cap, rna_sha256):
    cache = load_graph_cache(
        path, ("co", "sg"),
        {"cell_ids": cell_ids, "genes": list(genes), "manifest_sha": manifest_sha,
         "selection_seed": selection_seed, "pool_cap": pool_cap, "rna_sha256": rna_sha256},
        "PBMC scGPT", len(genes))
    return None if cache is None else (cache["co"], cache["sg"])


def write_scgpt_cache(path, co, sg, cell_ids, genes, manifest_sha, selection_seed, pool_cap,
                       rna_sha256):
    write_graph_cache(path, co=co, sg=sg, cell_ids=cell_ids, genes=genes,
                      manifest_sha=manifest_sha, selection_seed=selection_seed,
                      pool_cap=pool_cap, rna_sha256=rna_sha256)
