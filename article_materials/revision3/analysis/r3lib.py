"""Shared helpers for the Revision 3 analyses of ARRAY-D-26-05613.

The component randomization below is a pair-explicit copy of the registered
`screg_validation.batched_null.batched_components`. With the full non-self pair set,
restricted/unrestricted strata, seed 20261001 and own/none degree it must reproduce the
registered rows exactly (checked by `r3_reproduce_registered.py`). The copy adds three
things only: an explicit evaluated pair set, several degree designs evaluated on the same
draws, and score vectors that need not be model graphs (the co-expression comparator).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from r3paths import CORR, EXT, REV, VAL, W  # noqa: E402,F401
sys.path.insert(0, str(VAL))
from screg_validation.batched_null import basis, residual  # noqa: E402
from screg_validation.statistics import adjust, partial_spearman  # noqa: E402

READOUTS = ['geneformer_embed', 'random_init_geneformer', 'geneformer_attention',
            'scgpt_encoder', 'uce_encoder', 'scfoundation_encoder']
GENE_LEVEL = [r for r in READOUTS if r != 'geneformer_attention']
JOBS = ['brain_selected_panel-brain', 'brain_selected_panel-pbmc', 'historical_panel_corrected_source-brain']
LABEL = {'brain_selected_panel-brain': 'Brain panel / brain',
         'brain_selected_panel-pbmc': 'Brain panel / PBMC',
         'historical_panel_corrected_source-brain': 'Historical panel / brain'}
SEED = 20261001
KNOT_Q = (0.05, 0.275, 0.5, 0.725, 0.95)
DESIGNS = ['none', 'own', 'outdeg_only', 'indeg_only', 'log_degree', 'rank_degree',
           'spline_degree', 'strength', 'external_fibro', 'external_other']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def load_json(p):
    return json.loads(Path(p).read_text())


def write_json(p, obj):
    Path(p).write_text(json.dumps(obj, indent=1, default=float) + '\n')


def load_context(job, root=CORR, readouts=READOUTS):
    panels = load_json(root / 'panels.json')
    panel = panels[job.rsplit('-', 1)[0]]
    genes, tfg = panel['genes'], panel['tf_genes']
    a = dict(np.load(root / job / 'graphs.npz', allow_pickle=False))
    scores = {r: np.load(root / job / r / 'graph.npz', allow_pickle=False)['graph'] for r in readouts}
    return {'job': job, 'genes': list(genes), 'tf_genes': list(tfg),
            'tf': np.array([genes.index(t) for t in tfg]), 'a': a, 'scores': scores}


def missing_genes(graph):
    """Genes whose every off-diagonal score is undefined (gene-level missingness)."""
    g = np.array(graph, dtype=float, copy=True)
    np.fill_diagonal(g, np.nan)
    row = ~np.isfinite(g).any(1)
    col = ~np.isfinite(g).any(0)
    return row & col


def common_support(ctx, readouts=READOUTS, pair_readouts=None):
    """S: positive-variance genes with finite gene-level scores in every readout.
    E: non-self TF-target pairs within S finite for every readout in `pair_readouts`."""
    a, n = ctx['a'], len(ctx['genes'])
    keep = a['rna_variance'] > 0
    reasons = {'zero_rna_variance': [ctx['genes'][i] for i in np.flatnonzero(~keep)]}
    for r in readouts:
        if r == 'geneformer_attention':
            continue
        miss = missing_genes(ctx['scores'][r])
        reasons['missing_' + r] = [ctx['genes'][i] for i in np.flatnonzero(miss & keep)]
        keep &= ~miss
    S = np.flatnonzero(keep)
    pos = -np.ones(n, int)
    pos[S] = np.arange(len(S))
    tf_full = np.array([t for t in ctx['tf'] if keep[t]])
    tf_S = pos[tf_full]
    m = len(S)
    ii, jj = np.repeat(tf_S, m), np.tile(np.arange(m), len(tf_S))
    ok = ii != jj
    ii, jj = ii[ok], jj[ok]
    finite = np.ones(len(ii), bool)
    for r in (pair_readouts if pair_readouts is not None else readouts):
        sub = ctx['scores'][r][np.ix_(S, S)]
        finite &= np.isfinite(sub[ii, jj])
    co = a['coexpression'][np.ix_(S, S)]
    if not np.isfinite(co[ii, jj]).all():
        raise RuntimeError('comparator undefined inside common support')
    registered_pairs = len(ctx['tf']) * n - len(ctx['tf'])
    info = {'genes_panel': n, 'tfs_panel': int(len(ctx['tf'])), 'genes_S': int(m), 'tfs_S': int(len(tf_S)),
            'pairs_registered': int(registered_pairs), 'pairs_S': int(len(ii)), 'pairs_E': int(finite.sum()),
            'fraction_registered_pairs_retained': float(finite.sum() / registered_pairs),
            'excluded_genes': reasons}
    return S, tf_S, ii[finite], jj[finite], info


def subset(ctx, S):
    a = ctx['a']
    sub = {k: (a[k][np.ix_(S, S)] if a[k].ndim == 2 else a[k][S]) for k in a}
    scores = {r: ctx['scores'][r][np.ix_(S, S)] for r in ctx['scores']}
    return sub, scores


def strata_for(sub, tf_S, restricted):
    m = len(sub['peak_count'])
    if not restricted:
        return np.zeros(m, int)
    ob = np.digitize(sub['peak_count'], np.unique(np.quantile(sub['peak_count'], [.2, .4, .6, .8])))
    gb = np.digitize(sub['gc'], np.unique(np.quantile(sub['gc'], [1 / 3, 2 / 3])))
    return 15 * np.isin(np.arange(m), tf_S) + 3 * ob + gb


def ncs_basis(x, knots):
    """Natural cubic spline basis without intercept (as in the historical degree grid)."""
    K = len(knots)
    xi = np.asarray(knots, float)

    def d(k):
        return (np.clip(x - xi[k], 0, None) ** 3 - np.clip(x - xi[K - 1], 0, None) ** 3) / (xi[K - 1] - xi[k])
    return np.column_stack([x] + [d(k) - d(K - 2) for k in range(K - 2)])


def pair_ranks(values, index):
    """Average ranks of values[index] (ties averaged, 1-based), from gene-level values."""
    values = np.asarray(values)
    mult = np.bincount(index, minlength=len(values)).astype(float)
    present = np.flatnonzero(mult > 0)
    uv, inv = np.unique(values[present], return_inverse=True)
    tot = np.bincount(inv, weights=mult[present])
    before = np.concatenate([[0.0], np.cumsum(tot)[:-1]])
    rank_gene = np.zeros(len(values))
    rank_gene[present] = (before + (tot + 1) / 2)[inv]
    return rank_gene[index]


class Design:
    """Degree block for one design: fixed external columns and/or per-draw columns."""

    def __init__(self, name, ii, jj, graph, external=None):
        self.name, self.ii, self.jj = name, ii, jj
        self.fixed_extra = None
        self.knots = None
        if name in ('external_fibro', 'external_other', 'reference'):
            od, ind = external
            self.fixed_extra = np.column_stack([np.asarray(od, float)[ii], np.asarray(ind, float)[jj]])
        if name == 'spline_degree':
            od, ind = (graph > 0).sum(1), (graph > 0).sum(0)
            self.knots = (np.quantile(np.log1p(od[ii]), KNOT_Q), np.quantile(np.log1p(ind[jj]), KNOT_Q))

    def changing(self, g):
        """Per-draw degree block. Gene-level transforms are applied before expanding to pairs,
        which is elementwise identical to transforming the pair vectors (ranks use pair
        multiplicities, so they equal scipy average ranks over pairs exactly)."""
        ii, jj, n = self.ii, self.jj, self.name
        if n in ('none', 'external_fibro', 'external_other', 'reference'):
            return None
        od, ind = (g > 0).sum(1), (g > 0).sum(0)
        if n == 'own':
            return np.column_stack([od[ii], ind[jj]])
        if n == 'outdeg_only':
            return od[ii][:, None].astype(float)
        if n == 'indeg_only':
            return ind[jj][:, None].astype(float)
        if n == 'log_degree':
            return np.column_stack([np.log1p(od)[ii], np.log1p(ind)[jj]])
        if n == 'rank_degree':
            return np.column_stack([pair_ranks(od, ii), pair_ranks(ind, jj)])
        if n == 'spline_degree':
            return np.column_stack([ncs_basis(np.log1p(od.astype(float)), self.knots[0])[ii],
                                    ncs_basis(np.log1p(ind.astype(float)), self.knots[1])[jj]])
        if n == 'strength':
            return np.column_stack([g.sum(1)[ii], g.sum(0)[jj]])
        raise ValueError(n)


BASIS_IMPL = 'svd'


def qr_basis(values, *, reference_scale=0.):
    """Column-space basis identical in definition to `batched_null.basis`, computed as a reduced
    QR followed by an SVD of the small R factor (same singular values and rank tolerance)."""
    values = np.asarray(values, dtype=float)
    if values.shape[1] == 0:
        return np.empty((values.shape[0], 0))
    q, r = np.linalg.qr(values, mode='reduced')
    u, singular, _ = np.linalg.svd(r)
    if not len(singular) or singular[0] == 0:
        return np.empty((values.shape[0], 0))
    keep = singular > np.finfo(float).eps * max(values.shape) * max(singular[0], reference_scale)
    return q @ u[:, keep]


def spectral_norm(values):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(max(np.linalg.eigvalsh(values.T @ values).max(), 0.0)))


def components(score_vecs, graph, tf, ii, jj, controls, *, strata, permutations, seed, designs,
               external=None, backend='cpu', progress_every=1000, tag='', start=0, stop=None):
    """Shared M/D randomizations over an explicit pair set; one result per design x readout.

    RNG consumption matches `batched_components` exactly: per draw, one permutation of each
    stratum block (M), then for each TF row and block one permutation of its non-self targets (D).
    """
    graph = np.asarray(graph, float)
    n = len(graph)
    X = np.asarray(score_vecs, float)
    if X.ndim == 1:
        X = X[None, :]
    if X.shape[1] != len(ii) or not np.isfinite(X).all() or not np.isfinite(graph).all():
        raise RuntimeError('score or proxy not finite on evaluated pairs')
    controls = np.asarray(controls, float)
    blocks = [np.flatnonzero(strata == g) for g in np.unique(strata)]
    if not any(len(b) > 1 for b in blocks):
        raise RuntimeError('all strata singleton')
    base = np.column_stack([np.ones(len(ii)), controls])
    plans = []
    for name in designs:
        d = Design(name, ii, jj, graph, external=(external or {}).get(name))
        fixed = base if d.fixed_extra is None else np.column_stack([base, d.fixed_extra])
        q = basis(fixed)
        ranked = np.column_stack([rankdata(x) for x in X])
        if np.any(np.ptp(ranked, axis=0) == 0):
            raise RuntimeError('constant readout statistic')
        ranked -= ranked.mean(0, keepdims=True)
        xr = residual(ranked, q)
        norms = np.linalg.norm(xr, axis=0)
        tol = np.finfo(float).eps * max(fixed.shape) * 8
        if np.any(norms <= tol * np.linalg.norm(ranked, axis=0)):
            raise RuntimeError('fully conditioned readout')
        xr /= norms
        plans.append((d, q, xr, tol))

    if backend == 'cuda':
        import torch
        torch.backends.cuda.matmul.allow_tf32 = False
        dev = torch.device('cuda')
        tplans = [(d, torch.as_tensor(q, dtype=torch.float64, device=dev), torch.as_tensor(xr, dtype=torch.float64, device=dev), tol)
                  for d, q, xr, tol in plans]
        eps = float(np.finfo(float).eps)

        def vectors(g):
            raw = g[ii, jj]
            if np.ptp(raw) == 0:
                raise RuntimeError('constant randomized proxy')
            rp = rankdata(raw)
            rp -= rp.mean()
            rpt = torch.as_tensor(rp, dtype=torch.float64, device=dev)
            out = []
            for d, qt, xrt, tol in tplans:
                ch = d.changing(g)
                y = rpt - qt @ (qt.T @ rpt)
                if ch is None:
                    qc = torch.empty((len(ii), 0), dtype=torch.float64, device=dev)
                else:
                    cht = torch.as_tensor(np.asarray(ch, float), dtype=torch.float64, device=dev)
                    ref = float(torch.sqrt(torch.clamp(torch.linalg.eigvalsh(cht.T @ cht).max(), min=0)))
                    rc = cht - qt @ (qt.T @ cht)
                    qq, rr = torch.linalg.qr(rc, mode='reduced')
                    u, s, _ = torch.linalg.svd(rr)
                    s0 = float(s[0]) if len(s) else 0.0
                    if s0 == 0:
                        qc = torch.empty((len(ii), 0), dtype=torch.float64, device=dev)
                    else:
                        keep = s > eps * max(rc.shape) * max(s0, ref)
                        qc = qq @ u[:, keep]
                    y = y - qc @ (qc.T @ y)
                nrm = float(torch.linalg.vector_norm(y))
                if nrm <= tol * float(np.linalg.norm(rp)):
                    raise RuntimeError('fully conditioned randomized proxy')
                c = (xrt.T @ torch.column_stack([y / nrm, qc])).cpu().numpy()
                den = np.sqrt(np.maximum(0, 1 - np.sum(c[:, 1:] ** 2, axis=1)))
                if np.any(den < 1e-12):
                    raise RuntimeError('fully conditioned readout under randomization')
                out.append(c[:, 0] / den)
            return np.array(out)
    else:
        vectors = None

    def cpu_vectors(g):
        raw = g[ii, jj]
        if np.ptp(raw) == 0:
            raise RuntimeError('constant randomized proxy')
        rp = rankdata(raw)
        rp -= rp.mean()
        out = []
        for d, q, xr, tol in plans:
            ch = d.changing(g)
            if ch is None:
                qc = np.empty((len(ii), 0))
            elif BASIS_IMPL == 'svd':
                qc = basis(residual(ch, q), reference_scale=np.linalg.norm(ch, 2))
            else:
                qc = qr_basis(residual(ch, q), reference_scale=spectral_norm(ch))
            y = residual(residual(rp, q), qc)
            nrm = np.linalg.norm(y)
            if nrm <= tol * np.linalg.norm(rp):
                raise RuntimeError('fully conditioned randomized proxy')
            c = xr.T @ np.column_stack([y / nrm, qc])
            den = np.sqrt(np.maximum(0, 1 - np.sum(c[:, 1:] ** 2, axis=1)))
            if np.any(den < 1e-12):
                raise RuntimeError('fully conditioned readout under randomization')
            out.append(c[:, 0] / den)
        return np.array(out)  # designs x readouts

    if vectors is None:
        vectors = cpu_vectors
    observed = vectors(graph)
    stop = permutations if stop is None else stop
    nulls = {'M': np.empty((stop - start,) + observed.shape), 'D': np.empty((stop - start,) + observed.shape)}
    rng = np.random.default_rng(seed)
    t0 = time.time()
    for k in range(stop):
        perm = np.arange(n)
        for b in blocks:
            perm[b] = rng.permutation(b)
        if k >= start:
            nulls['M'][k - start] = vectors(graph[np.ix_(perm, perm)])
            shuffled = graph.copy()
            for t in tf:
                for b in blocks:
                    el = b[b != t]
                    shuffled[t, el] = graph[t, rng.permutation(el)]
            nulls['D'][k - start] = vectors(shuffled)
        else:
            for t in tf:
                for b in blocks:
                    rng.permutation(b[b != t])
        if progress_every and k >= start and (k + 1) % progress_every == 0:
            log(tag, 'draw', k + 1, '/', stop, f'{(time.time() - t0) / (k + 1 - start):.3f}s/draw')
    if start or stop != permutations:
        return {'observed': observed, 'nulls': nulls, 'start': start, 'stop': stop, 'stratum_sizes': [len(b) for b in blocks]}
    pM = (1 + (np.abs(nulls['M']) >= np.abs(observed)[None]).sum(0)) / (permutations + 1)
    pD = (1 + (np.abs(nulls['D']) >= np.abs(observed)[None]).sum(0)) / (permutations + 1)
    return {'observed': observed, 'pM': pM, 'pD': pD, 'pIUT': np.maximum(pM, pD), 'nulls': nulls,
            'stratum_sizes': [len(b) for b in blocks]}


def tf_cluster_bootstrap(xvec, yvec, cvec, controls, ii, resamples=2000, seed=20261031):
    """Percentile interval for rho(x,y|C) - rho(c,y|C), resampling TF rows."""
    rows = np.unique(ii)
    idx = {t: np.flatnonzero(ii == t) for t in rows}
    rng = np.random.default_rng(seed)
    point = partial_spearman(xvec, yvec, controls) - partial_spearman(cvec, yvec, controls)
    boots = np.empty(resamples)
    for b in range(resamples):
        pick = rng.choice(rows, size=len(rows), replace=True)
        sel = np.concatenate([idx[t] for t in pick])
        boots[b] = (partial_spearman(xvec[sel], yvec[sel], controls[sel])
                    - partial_spearman(cvec[sel], yvec[sel], controls[sel]))
    lo, hi = np.quantile(boots, [.025, .975])
    return float(point), [float(lo), float(hi)]
