"""Preregistered structural null stress, separate from biological validation."""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .common import ContractError
from .batched_null import batched_components


def wilson(k, n):
    z = float(norm.ppf(.975))
    p, d = k / n, 1 + z * z / n
    center = (p + z * z / (2 * n)) / d
    radius = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(max(0., center-radius)), float(min(1., center+radius))]


def stress(graph, tf_indices, controls, strata, *, replicates, permutations, seed, alpha=.05,
           backend="cpu", gpu_policy=None, replicate_batch=8, check=None):
    if replicates < 1 or not 0 < alpha < 1:
        raise ContractError("invalid null stress registration")
    g = np.asarray(graph, float)
    n = len(g)
    groups = np.asarray(strata)
    if groups.shape != (n,):
        raise ContractError("stress strata shape mismatch")
    rng = np.random.default_rng(seed)
    ii, jj = np.repeat(tf_indices, n), np.tile(np.arange(n), len(tf_indices))
    ii, jj = ii[ii != jj], jj[ii != jj]
    rows = []
    for scenario in ("iid_edge", "stratum_heteroskedastic_shared_shock", "stratum_mean_misspecification"):
        counts = {"M": 0, "D": 0, "both": 0}
        mean_rho = []
        for start in range(0, replicates, replicate_batch):
            if check:
                check()
            scores = []
            for _ in range(start, min(start + replicate_batch, replicates)):
                noise = rng.normal(size=(n, n))
                if scenario != "iid_edge":
                    scales = 1 + (groups - groups.min()) / max(1, groups.max() - groups.min())
                    noise += rng.normal(size=n)[:, None] * scales[:, None]
                    noise += rng.normal(size=n)[None, :] * scales[None, :]
                    if scenario == "stratum_mean_misspecification":
                        noise += scales[:, None] + scales[None, :]
                scores.append(noise)
            tested = batched_components(scores, g, tf_indices, controls, strata=groups,
                        permutations=permutations, seed=seed, degree_mode="own",
                        backend=backend, gpu_policy=gpu_policy, check=check)
            for result in tested:
                counts["M"] += result["pM"] <= alpha
                counts["D"] += result["pD"] <= alpha
                counts["both"] += max(result["pM"], result["pD"]) <= alpha
                mean_rho.append(result["observed"])
        rows.append({"scenario": scenario, "replicates": replicates,
                     "mean_observed": float(np.mean(mean_rho)),
                     "is_type_I_calibration": scenario != "stratum_mean_misspecification",
                     "rates": {k: {"rejections": int(v), "proportion": v / replicates,
                                    "wilson_ci95": wilson(v, replicates)} for k, v in counts.items()}})
    return {"rows": rows, "alpha": alpha, "permutations": permutations, "seed": seed,
            "randomization_stream": "shared fixed stream for all replicates; intervals conditional on this stream",
            "scope": "registered synthetic mechanisms only; shared-shock null validity additionally requires exchangeability under these strata; no unconditional real-data certification"}
