"""Shared randomizations and small BLAS batches across readouts; no giant null bank."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from scipy.stats import rankdata

from .common import ContractError, digest, load_json, write_json
from .hardware import require_gpu


def basis(values, *, reference_scale=0.):
    values = np.asarray(values, dtype=float)
    u, singular, _ = np.linalg.svd(values, full_matrices=False)
    if not len(singular) or singular[0] == 0:
        return np.empty((values.shape[0], 0))
    return u[:, singular > np.finfo(float).eps * max(values.shape) * max(singular[0], reference_scale)]


def residual(values, q):
    return values - q @ (q.T @ values)


def batched_components(scores, graph, tf_indices, controls, *, strata, permutations, seed,
                       degree_mode, reference_degree=None, backend="cpu", gpu_policy=None,
                       batch_permutations=8, max_batch_bytes=256 * 1024 ** 2, check=None,
                       checkpoint=None, checkpoint_context=None, checkpoint_every=128):
    """FWL equivalent to direct partial Spearman; changing own-degree refits per draw.

    Gene/row randomizations are shared across all readouts. Matrix multiplication
    batches have a fixed memory cap; output stores only per-readout null statistics.
    """
    graph, scores = np.asarray(graph, float), np.asarray(scores, float)
    n = len(graph)
    if scores.ndim == 2:
        scores = scores[None, :, :]
    if graph.shape != (n, n) or scores.shape[1:] != graph.shape or not np.isfinite(graph).all():
        raise ContractError("batched randomization graph shape")
    tf = np.asarray(tf_indices, int)
    if not len(tf) or len(set(tf)) != len(tf) or tf.min() < 0 or tf.max() >= n:
        raise ContractError("invalid TF rows")
    ii, jj = np.repeat(tf, n), np.tile(np.arange(n), len(tf))
    ii, jj = ii[ii != jj], jj[ii != jj]
    if not np.isfinite(scores[:, ii, jj]).all():
        raise ContractError("unobserved evaluated score; frozen pairs cannot be zero imputed")
    controls = np.asarray(controls, float)
    if controls.ndim != 2 or controls.shape[0] != len(ii) or not np.isfinite(controls).all():
        raise ContractError("batched nuisance shape")
    groups = np.asarray(strata)
    if groups.shape != (n,) or permutations < 1 or degree_mode not in {"own", "reference", "none"}:
        raise ContractError("invalid randomization registration")
    blocks = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    if not any(len(b) > 1 for b in blocks):
        raise ContractError("all strata singleton")
    fixed = np.column_stack([np.ones(len(ii)), controls])
    if degree_mode == "reference":
        if reference_degree is None or any(np.asarray(v).shape != (n,) for v in reference_degree):
            raise ContractError("reference degree shape")
        fixed = np.column_stack([fixed, reference_degree[0][ii], reference_degree[1][jj]])
    q = basis(fixed)
    ranked = np.column_stack([rankdata(s[ii, jj]) for s in scores])
    if np.any(np.ptp(ranked, axis=0) == 0):
        raise ContractError("constant readout statistic")
    ranked -= ranked.mean(0, keepdims=True)
    xr = residual(ranked, q)
    norms = np.linalg.norm(xr, axis=0)
    tolerance = np.finfo(float).eps * max(fixed.shape) * 8
    if np.any(norms <= tolerance * np.linalg.norm(ranked, axis=0)):
        raise ContractError("constant or fully conditioned readout")
    xr /= norms
    def vector(g):
        if degree_mode == "own":
            changing = np.column_stack([(g > 0).sum(1)[ii], (g > 0).sum(0)[jj]])
            # Constant degree columns project to floating-point noise. Their rank
            # tolerance must use the original column scale, not that tiny residue.
            qc = basis(residual(changing, q), reference_scale=np.linalg.norm(changing, 2))
        else:
            qc = np.empty((len(ii), 0))
        raw = g[ii, jj]
        if np.ptp(raw) == 0:
            raise ContractError("constant or fully conditioned randomized proxy")
        ranked_proxy = rankdata(raw)
        ranked_proxy -= ranked_proxy.mean()
        y = residual(residual(ranked_proxy, q), qc)
        norm = np.linalg.norm(y)
        if norm <= tolerance * np.linalg.norm(ranked_proxy):
            raise ContractError("constant or fully conditioned randomized proxy")
        return np.column_stack([y / norm, qc])
    observed_vector = vector(graph)
    cross = xr.T @ observed_vector
    den = np.sqrt(np.maximum(0, 1 - np.sum(cross[:, 1:] ** 2, axis=1)))
    if np.any(den < 1e-12):
        raise ContractError("fully conditioned observed readout")
    observed = cross[:, 0] / den
    if backend not in {"cpu", "cuda"}:
        raise ContractError("explicit null backend required")
    if backend == "cuda":
        import torch
        require_gpu(gpu_policy or {})
        torch.cuda.set_device(int((gpu_policy or {}).get("index", 0)))
        torch.backends.cuda.matmul.allow_tf32 = False
        device_scores = torch.as_tensor(xr.T, dtype=torch.float64, device="cuda")
    max_columns = 3 if degree_mode == "own" else 1
    batch_size = max(1, min(batch_permutations, max_batch_bytes // (len(ii) * max_columns * 8)))
    nulls = {"M": np.empty((permutations, len(scores))), "D": np.empty((permutations, len(scores)))}
    rng = np.random.default_rng(seed)
    next_draw = 0
    if checkpoint is not None:
        checkpoint = Path(checkpoint)
        if checkpoint_context is None:
            raise ContractError("draw checkpoint requires an explicit stage context")
        if checkpoint.exists():
            state = load_json(checkpoint)
            contents = {k: v for k, v in state.items() if k != "state_sha256"}
            if state.get("state_sha256") != digest(contents) or state.get("context_sha256") != digest(checkpoint_context):
                raise ContractError("draw checkpoint content or context drift")
            next_draw = state["next_draw"]
            if type(next_draw) is not int or not 0 <= next_draw <= permutations:
                raise ContractError("invalid draw checkpoint position")
            for kind in nulls:
                saved = np.asarray(state["nulls"][kind], float).reshape(next_draw, len(scores))
                if not np.isfinite(saved).all():
                    raise ContractError("nonfinite saved randomization statistics")
                nulls[kind][:next_draw] = saved
            rng.bit_generator.state = state["rng_state"]
    saved_through = next_draw
    pending = []
    def flush():
        nonlocal saved_through
        if check:
            check()
        matrix = np.column_stack([v for _, _, v in pending])
        if backend == "cuda":
            require_gpu(gpu_policy or {})
            packed = torch.as_tensor(matrix, device="cuda", dtype=torch.float64)
            values = (device_scores @ packed).cpu().numpy()
            del packed
        else:
            values = xr.T @ matrix
        offset = 0
        for kind, index, v in pending:
            block = values[:, offset:offset + v.shape[1]]
            denominator = np.sqrt(np.maximum(0, 1 - np.sum(block[:, 1:] ** 2, axis=1)))
            if np.any(denominator < 1e-12):
                raise ContractError("fully conditioned randomized readout")
            nulls[kind][index] = block[:, 0] / denominator
            offset += v.shape[1]
        pending.clear()
        if checkpoint is not None and (k + 1 - saved_through >= checkpoint_every or k + 1 == permutations):
            state = {"context_sha256": digest(checkpoint_context), "next_draw": k + 1,
                     "rng_state": rng.bit_generator.state,
                     "nulls": {name: values[:k + 1].tolist() for name, values in nulls.items()}}
            write_json(checkpoint, {**state, "state_sha256": digest(state)})
            saved_through = k + 1
    for k in range(next_draw, permutations):
        permutation = np.arange(n)
        for block in blocks:
            permutation[block] = rng.permutation(block)
        pending.append(("M", k, vector(graph[np.ix_(permutation, permutation)])))
        shuffled = graph.copy()
        for t in tf:
            for block in blocks:
                eligible = block[block != t]
                shuffled[t, eligible] = graph[t, rng.permutation(eligible)]
        pending.append(("D", k, vector(shuffled)))
        if len(pending) >= batch_size or k == permutations - 1:
            flush()
    if backend == "cuda":
        del device_scores
        torch.cuda.empty_cache()
    results = []
    for index in range(len(scores)):
        pm, pd = [(1 + (np.abs(nulls[k][:, index]) >= abs(observed[index])).sum()) / (permutations + 1)
                  for k in ("M", "D")]
        results.append({"observed": float(observed[index]), "pM": float(pm), "pD": float(pd), "pIUT": float(max(pm, pd)),
                        "permutations": permutations, "degree_mode": degree_mode,
                        "stratum_sizes": [len(b) for b in blocks], "backend": backend,
                        "nulls": {k: nulls[k][:, index].tolist() for k in nulls},
                        "scope": "shared conditional component reference distributions; strata exchangeability remains an assumption"})
    return results
