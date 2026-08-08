#!/usr/bin/env python
"""Shared co-expression baseline runner, parameterised by tissue.

The baseline is scored at its own primary rung (no self-partial, use_coexp=False,
full confounds) because partialling co-expression out of co-expression is
degenerate by construction. Same edge set, confounds, and null machinery as the
pooled audit rows. Explicit integer seeds make both Monte Carlo tests
reproducible across Python processes.
"""
import os
import time
from typing import Optional

import numpy as np

import audit_utils as au
import fixed_panel_audit as fpa


N_PERM = 999
# The four gene-identity confounds (peak count, gene length, detection, GC) are
# properties of gene identity, so both tissues read the same cached vectors.
CONFOUND_CACHE = "pbmc_confounds_v2.npz"


def run_baseline(tissue: str, atac_file: str, proxy_npz: str, fm_npz: str,
                 mantel_seed: int, degree_seed: int, out_name: str,
                 n_perm: int = N_PERM,
                 confound_cache: Optional[str] = None) -> dict:
    """Score the co-expression baseline for one tissue through both nulls.

    Writes ``results/v2/<out_name>`` and returns the document that was written.
    """
    G, tf, _types = fpa.load_consensus_proxy(f"{fpa.OUT}/{proxy_npz}")
    co = np.load(f"{fpa.OUT}/{fm_npz}")["co"].astype(np.float32)
    cache = np.load(f"{fpa.OUT}/{confound_cache or CONFOUND_CACHE}", allow_pickle=False)
    gl, dv, gc = cache["genelen"], cache["detv"], cache["gc"]
    au.log(f"{tissue} peakcount...")
    pc = fpa.peak_counts(atac_file)
    od, ind = fpa.graph_degrees(G)
    ii, jj = fpa.panel_edge_indices(tissue, tf, G.shape[0])
    co_v, at_v = co[ii, jj], G[ii, jj]
    au.log(f"{tissue} edges {len(ii)}")
    obs = fpa.partial_rho_obs_sliced(co_v, at_v, co_v, jj, ii, pc, gl, dv, gc,
                                     od, ind, False, "full")
    au.log(f"{tissue} coexp baseline obs (degree_only, no self-partial): {obs:.6f}")
    t = time.time()
    man = fpa.mantel_randomization(co_v, at_v, co_v, jj, ii, pc, gl, dv, gc,
                                   od, ind, G, False, "full", obs, n_perm,
                                   seed=mantel_seed)
    deg = fpa.degree_preserving_null(co_v, at_v, co_v, jj, ii, pc, gl, dv, gc,
                                     od, ind, G, np.unique(tf), False, "full", obs, n_perm,
                                     seed=degree_seed)
    au.log(f"{tissue} baseline: rho={obs:+.6f} pM={man['p_mc']} pD={deg['p_mc']} "
           f"zM={man['z']:.2f} zD={deg['z']:.2f}  ({time.time() - t:.0f}s)")
    document = {
        "schema_version": 2,
        "tissue": tissue, "model_label": "co_expression_baseline",
        "rung": "degree_only_no_selfpartial", "observed_rho": float(obs),
        "pM": man["p_mc"], "pD": deg["p_mc"], "zM": man["z"], "zD": deg["z"],
        "n_perm": n_perm, "mantel_seed": mantel_seed,
        "degree_seed": degree_seed, "seed_contract": "explicit_integer_v1",
    }
    au.write_json_atomic(f"{fpa.OUT}/{out_name}", document, indent=1, sort_keys=False,
                         allow_nan=True, newline=False)
    au.log(f"wrote {out_name}")
    return document


def data_root() -> str:
    """External data lake root used to locate ATAC inputs."""
    return os.environ.get(
        "SCREG_DATA_ROOT",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data"))
