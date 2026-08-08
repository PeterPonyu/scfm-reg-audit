#!/usr/bin/env python
"""Shared infrastructure helpers used by every audit script.

Collects the utilities that were previously copy-pasted across the pipeline:
timestamped logging, chunked file hashing, crash-safe (atomic) JSON/NPZ writes,
Benjamini-Hochberg q-values, and the plus-one Monte-Carlo p-value summary.
"""
import hashlib
import json
import os
import tempfile
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np


# ----------------------------- logging ---------------------------------------
def log(*values, flush: bool = True) -> None:
    """Print values prefixed with a wall-clock timestamp."""
    print(f"[{time.strftime('%H:%M:%S')}]", *values, flush=flush)


# ----------------------------- hashing ---------------------------------------
def sha256_file(path, chunk_size: int = 1024 * 1024) -> str:
    """SHA-256 of a file, read in chunks so large caches never load into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    """SHA-256 of an array's contiguous buffer."""
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# ----------------------------- atomic writes ---------------------------------
def atomic_temp_path(path, prefix: Optional[str] = None, suffix: str = ".tmp") -> str:
    """Create an empty temporary file next to ``path`` and return its name.

    Staging in the destination directory keeps the final ``os.replace`` atomic.
    """
    destination = Path(path)
    parent = destination.parent if str(destination.parent) else Path(".")
    parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=prefix if prefix is not None else destination.name,
        suffix=suffix, dir=parent)
    os.close(fd)
    return temporary


def write_json_atomic(path, document, indent: int = 2, sort_keys: bool = True,
                      allow_nan: bool = False, newline: bool = True) -> None:
    """Write JSON to a sibling temporary file, then rename it over ``path``."""
    temporary = atomic_temp_path(path, suffix=".json")
    try:
        with open(temporary, "w") as fh:
            json.dump(document, fh, indent=indent, sort_keys=sort_keys, allow_nan=allow_nan)
            if newline:
                fh.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_npz_atomic(path, **arrays) -> None:
    """Write an .npz to a sibling temporary file, then rename it over ``path``."""
    temporary = atomic_temp_path(path, suffix=".npz")
    try:
        np.savez(temporary, **arrays)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


# ----------------------------- BH-FDR ----------------------------------------
def bh_qvalues(pvalues: Iterable[float], round_to: Optional[int] = None) -> List[float]:
    """Benjamini-Hochberg q-values in the input order, clipped to [0, 1]."""
    p = np.asarray(list(pvalues), dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    monotone = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(monotone, 0, 1)
    if round_to is None:
        return [float(x) for x in out]
    return [round(float(x), round_to) for x in out]


# ----------------------------- Monte-Carlo p-values --------------------------
def mc_pvalue(null: np.ndarray, observed: float, n_perm: int) -> tuple:
    """Plus-one two-sided Monte-Carlo p-value.

    p_mc = (count(|null| >= |observed|) + 1) / (n_perm + 1). Returns (p_mc, count).
    """
    count = int(np.sum(np.abs(null) >= abs(observed)))
    return (count + 1) / (n_perm + 1), count


def mc_null_summary(null: np.ndarray, observed: float, n_perm: int,
                    seed: Optional[int]) -> Dict[str, float]:
    """Point estimate, plus-one p_mc, and null moments for one null distribution."""
    p_mc, count = mc_pvalue(null, observed, n_perm)
    return {
        "p_mc": float(p_mc),
        "N_perm": int(n_perm),
        "seed": int(seed) if seed is not None else None,
        "resolution": float(1 / (n_perm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "z": float((observed - null.mean()) / (null.std() + 1e-9)),
        "null_obs_count_at_or_above_obs": count,
    }


# ----------------------------- seeds -----------------------------------------
def spawn_int_seeds(seed_sequence: np.random.SeedSequence, n: int) -> List[int]:
    """Spawn n SeedSequence children and return n Python int seeds."""
    return [int(child.generate_state(1, dtype=np.uint64)[0])
            for child in seed_sequence.spawn(n)]
