#!/usr/bin/env python3
"""R3-1: calibration of unrestricted and restricted component tests under synthetic stress.

Extends the published component stress test (200 replicates per scenario, 999 draws, raw
p <= 0.05) without changing it. For the two historical contexts the published synthetic
readouts are regenerated from their recorded seeds and the unrestricted reference
distributions are regenerated from their recorded seeds; both must reproduce the saved
statistics before anything new is read. The restricted scheme (TF status x opportunity
quintile x GC tertile) is then applied to the same synthetic readouts, and a further
scenario scales the shared shocks by RNA detection, which is not a stratifier.
The genuine-brain context (brain-selected panel, GSE174367) uses the corrective covariates.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from r3paths import CORR, EXT, HIST_INPUTS, HIST_STRESS, REV  # noqa: E402,F401
sys.path.insert(0, str(EXT / 'analysis'))
from extend_injection_power import basis, null_vectors, residual, zscore  # noqa: E402

PUBLISHED_SEED = [20260927, 5613, 51]
NEW_SEED = [20261031, 5613]
N_PERM = 999
N_REP = 200


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def wilson(k, n, z=1.959963984540054):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


def nominal(k, n):
    lo, hi = wilson(k, n)
    rate = k / n
    return bool(lo <= 0.05 <= hi and 0.025 <= rate <= 0.075)


def load_context(name):
    """Return the arrays the published stress test used, plus stratifiers."""
    if name in ('brain', 'pbmc'):
        g = dict(np.load(HIST_INPUTS / f'{name}.npz', allow_pickle=False))
        seq = dict(np.load(HIST_INPUTS / 'sequence_opportunity.npz', allow_pickle=False))
        g['tissue'] = name
        ng = len(g['genes'])
        tf = g['tf']
        # Published shocks: TF and target sequence opportunity, sequence GC.
        a_tf = seq['tf_opportunity_outdeg']
        a_tg = seq['target_opportunity_indeg']
        a_noise = seq['sequence_gc']
        # Published restricted strata (restricted_randomization.py).
        ob = np.digitize(seq['target_opportunity_indeg'], np.unique(np.quantile(seq['target_opportunity_indeg'], [.2, .4, .6, .8])))
        gb = np.digitize(seq['sequence_gc'], np.unique(np.quantile(seq['sequence_gc'], [1 / 3, 2 / 3])))
        groups = ob * 3 + gb + 15 * np.isin(np.arange(ng), tf)
        excluded = np.setdiff1d(np.arange(ng), np.unique(np.r_[g['ii'], g['jj']]))
        groups[excluded] = 30 + np.arange(len(excluded))
        g['stratifier_note'] = 'TF status x target sequence-opportunity quintile x sequence-GC tertile; excluded markers singleton'
    elif name in EXTRA:
        return load_extra(name)
    else:
        panels = json.loads((CORR / 'panels.json').read_text())
        panel = panels['brain_selected_panel']
        genes, tfg = panel['genes'], panel['tf_genes']
        a = dict(np.load(CORR / 'brain_selected_panel-brain' / 'graphs.npz', allow_pickle=False))
        ng = len(genes)
        tf = np.array([genes.index(t) for t in tfg])
        ii, jj = np.repeat(tf, ng), np.tile(np.arange(ng), len(tf))
        keep = ii != jj
        ii, jj = ii[keep], jj[keep]
        if not np.isfinite(a['coexpression'][ii, jj]).all():
            raise SystemExit('genuine-brain comparator undefined on some pairs')
        g = {'G': a['proxy'].astype(np.float64), 'co': a['coexpression'][ii, jj], 'ii': ii, 'jj': jj, 'tf': tf,
             'pc': a['peak_count'].astype(np.float64), 'gl': a['gene_length'].astype(np.float64),
             'gc': a['gc'].astype(np.float64), 'det': a['detection'].astype(np.float64),
             'genes': np.array(genes), 'tissue': 'genuine_brain'}
        # Corrective covariates: shocks scaled by peak count and GC, the corrective stratifiers.
        a_tf = g['pc'][tf]
        a_tg = g['pc']
        a_noise = g['gc']
        ob = np.digitize(g['pc'], np.unique(np.quantile(g['pc'], [.2, .4, .6, .8])))
        gb = np.digitize(g['gc'], np.unique(np.quantile(g['gc'], [1 / 3, 2 / 3])))
        groups = 15 * np.isin(np.arange(ng), tf) + 3 * ob + gb
        g['stratifier_note'] = 'TF status x peak-count quintile x linked-window GC tertile (corrective restricted scheme)'
    g['groups'] = groups
    g['shock_A'] = (a_tf, a_tg, a_noise)
    g['shock_B'] = (g['det'][g['tf']], g['det'], g['det'])
    return g


EXTRA = {'new_pbmc': ('corr', 'brain_selected_panel-pbmc'), 'old_brain': ('corr', 'historical_panel_corrected_source-brain'),
         'bmmc_13272': ('bmmc', '13272'), 'bmmc_28483': ('bmmc', '28483')}


def load_extra(name):
    """Calibration contexts on the same edges and covariates as the reported analyses."""
    kind, key = EXTRA[name]
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import r3lib as RL
    if kind == 'corr':
        ctx = RL.load_context(key)
        S, tf, ii, jj, info = RL.common_support(ctx)
        sub, _ = RL.subset(ctx, S)
        g = {'G': sub['proxy'].astype(np.float64), 'co': sub['coexpression'][ii, jj], 'ii': ii, 'jj': jj, 'tf': tf,
             'pc': sub['peak_count'].astype(np.float64), 'gl': sub['gene_length'].astype(np.float64),
             'gc': sub['gc'].astype(np.float64), 'det': sub['detection'].astype(np.float64),
             'genes': np.array([ctx['genes'][i] for i in S]), 'tissue': name}
        g['groups'] = RL.strata_for(sub, tf, True)
        g['shock_A'] = (g['pc'][tf], g['pc'], g['gc'])
        g['stratifier_note'] = 'TF status x peak-count quintile x linked-window GC tertile on common-support edges'
    else:
        import r3_bmmc as BM
        plan, panel, prepared, conf = BM.setup()
        S, _detail = BM.common_genes(plan, panel, prepared)
        pos = -np.ones(len(panel['genes']), int)
        pos[S] = np.arange(len(S))
        tf = np.array([pos[panel['genes'].index(t)] for t in panel['tf_genes'] if pos[panel['genes'].index(t)] >= 0])
        ii, jj = BM.pairs_for(tf, len(S))
        a, var, scores = BM.donor_arrays(prepared, key)
        g = {'G': a['proxy'][np.ix_(S, S)].astype(np.float64), 'co': a['coexpression'][np.ix_(S, S)][ii, jj], 'ii': ii, 'jj': jj,
             'tf': tf, 'pc': conf['peak_count'][S].astype(np.float64), 'gl': conf['gene_length'][S].astype(np.float64),
             'gc': conf['gc'][S].astype(np.float64), 'det': conf['detection'][S].astype(np.float64),
             'genes': np.array([panel['genes'][i] for i in S]), 'tissue': name}
        g['groups'] = BM.restricted_strata(conf, S, tf)
        opp = conf['opportunity_degree'][S].astype(np.float64)
        g['shock_A'] = (opp[tf], opp, g['gc'])
        g['stratifier_note'] = 'TF status x opportunity-degree quintile x GC tertile (registered bone-marrow rule) on the common gene set'
    g['shock_B'] = (g['det'][g['tf']], g['det'], g['det'])
    return g


def design(g):
    ii, jj, G = g['ii'], g['jj'], g['G']
    od = (G > 0).sum(1).astype(np.float32)
    ind = (G > 0).sum(0).astype(np.float32)
    fixed = np.column_stack([np.ones(len(ii)), rankdata(g['co']), *[zscore(g[k][jj]) for k in ('pc', 'gl', 'det', 'gc')]])
    X = np.column_stack([fixed, zscore(od[ii]), zscore(ind[jj])])
    q = basis(X)
    y = residual(rankdata(G[ii, jj]), q)
    y -= y.mean()
    y /= np.linalg.norm(y)
    return fixed, od, ind, q, y


def synthetic(g, kind, seeds, a):
    """Synthetic readouts exactly as the published stress test builds them."""
    ii, jj = g['ii'], g['jj']
    n = len(ii)
    tf = g['tf']
    ntf, ng = len(tf), len(g['genes'])
    rowi = np.searchsorted(tf, ii)
    if kind != 'iid':
        a_tf, a_tg, a_noise = a
        s_tf = 1 + 2 * rankdata(a_tf) / (ntf + 1)
        s_tg = 1 + 2 * rankdata(a_tg) / (ng + 1)
        scale = 1 + abs(zscore(a_noise))
    out = np.empty((len(seeds), n))
    for k, seed in enumerate(seeds):
        rng = np.random.default_rng(int(seed))
        e = rng.normal(size=n)
        if kind == 'iid':
            v = e
        elif kind == 'mean':
            raise ValueError('mean stress handled separately')
        else:
            v = (s_tf[rowi] * rng.normal(size=ntf)[rowi] + s_tg[jj] * rng.normal(size=ng)[jj] + scale[jj] * e) / np.sqrt(3)
        out[k] = rankdata(v)
    return out


def mean_stress(g, seeds):
    seq = dict(np.load(HIST_INPUTS / 'sequence_opportunity.npz', allow_pickle=False))
    ii, jj = g['ii'], g['jj']
    rowi = np.searchsorted(g['tf'], ii)
    opp = seq['opportunity_counts'][rowi, jj]
    mean = zscore(np.log1p(opp))
    out = np.empty((len(seeds), len(ii)))
    for k, seed in enumerate(seeds):
        rng = np.random.default_rng(int(seed))
        out[k] = rankdata(mean + rng.normal(size=len(ii)))
    return out


def observed_stats(x, q, y):
    rx = x - (x @ q) @ q.T
    rx -= rx.mean(1, keepdims=True)
    return (rx @ y) / np.linalg.norm(rx, axis=1)


def restricted_vectors(g, fixed, od, ind, nt, seed):
    """Restricted analogue of the published null_vectors generator."""
    ii, jj, G = g['ii'], g['jj'], g['G']
    ng = G.shape[0]
    groups = g['groups']
    blocks = [np.flatnonzero(groups == k) for k in np.unique(groups)]
    xx = fixed if nt == 'mantel' else np.column_stack([fixed, zscore(od[ii])])
    q = basis(xx)
    tf = np.unique(g['tf'])
    colblocks = {int(t): [b[b != t] for b in blocks if np.any(b != t)] for t in tf}
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    for _ in range(N_PERM):
        if nt == 'mantel':
            p = np.arange(ng)
            for b in blocks:
                p[b] = rng.permutation(b)
            yv = rankdata(G[p[ii], p[jj]])
            ch = np.column_stack([zscore(od[p][ii]), zscore(ind[p][jj])])
        else:
            v = G.copy()
            for t in tf:
                for b in colblocks[int(t)]:
                    v[t, b] = G[t, rng.permutation(b)]
            yv = rankdata(v[ii, jj])
            newind = (v > 0).sum(0).astype(np.float64)
            ch = zscore(newind[jj])[:, None]
        qc = basis(residual(ch, q))
        yr = residual(residual(yv, q), qc)
        yr -= yr.mean()
        yield np.column_stack([yr / np.linalg.norm(yr), qc])


def project(xs, fixed, od, ii, nt, vectors, device):
    xx = fixed if nt == 'mantel' else np.column_stack([fixed, zscore(od[ii])])
    q = torch.as_tensor(basis(xx).T, device=device, dtype=torch.float64)
    xr = torch.as_tensor(xs, device=device, dtype=torch.float64)
    xr.sub_((xr @ q.T) @ q)
    xr.sub_(xr.mean(dim=1, keepdim=True))
    xr.div_(torch.linalg.vector_norm(xr, dim=1, keepdim=True))
    out = np.empty((N_PERM, xs.shape[0]))
    for k, v in enumerate(vectors):
        c = (xr @ torch.as_tensor(v, device=device, dtype=torch.float64)).cpu().numpy()
        den = np.sqrt(np.maximum(0, 1 - np.sum(c[:, 1:] ** 2, axis=1)))
        if np.any(den <= 0):
            raise RuntimeError('fully conditioned synthetic readout')
        out[k] = c[:, 0] / den
        if (k + 1) % 200 == 0:
            log(nt, k + 1, '/', N_PERM)
    del xr, q
    if device == 'cuda':
        torch.cuda.empty_cache()
    return out


def pvalues(nulls, obs):
    return (1 + np.sum(np.abs(nulls) >= np.abs(obs)[None, :], axis=0)) / (N_PERM + 1)


def summarize(pm, pd, labels, scenarios):
    rows = []
    for sc in scenarios:
        sel = labels == sc
        row = {'scenario': sc, 'n': int(sel.sum())}
        for lab, p in (('M', pm[sel]), ('D', pd[sel]), ('joint', np.maximum(pm[sel], pd[sel]))):
            k = int((p <= .05).sum())
            row[lab] = {'count': k, 'rate': k / int(sel.sum()), 'wilson95': wilson(k, int(sel.sum())),
                        'nominal_bradley_and_wilson': nominal(k, int(sel.sum()))}
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--context', choices=['brain', 'pbmc', 'genuine_brain', 'new_pbmc', 'old_brain', 'bmmc_13272', 'bmmc_28483'], required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--device', default='cuda')
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    g = load_context(a.context)
    tc = {'brain': 0, 'pbmc': 1, 'genuine_brain': 2, 'new_pbmc': 3, 'old_brain': 4, 'bmmc_13272': 5, 'bmmc_28483': 6}[a.context]
    fixed, od, ind, q, y = design(g)
    report = {'context': a.context, 'stratifiers': g['stratifier_note'], 'n_pairs': int(len(g['ii'])),
              'n_replicates_per_scenario': N_REP, 'n_permutations': N_PERM, 'threshold': 'raw p <= 0.05',
              'nominal_rule': 'Wilson 95% interval contains 0.05 and rate in [0.025, 0.075] (Bradley liberal criterion)',
              'plan_sha256': (Path(__file__).parent / 'R3_PLAN.sha256').read_text().split()[0]}
    if a.context in ('brain', 'pbmc'):
        seeds = [int(s.generate_state(1, dtype=np.uint64)[0]) for s in np.random.SeedSequence(PUBLISHED_SEED + [tc]).spawn(600)]
        saved = dict(np.load(HIST_STRESS / a.context / 'results.npz', allow_pickle=False))
        if not np.array_equal(np.array(seeds, dtype=np.uint64), saved['seed']):
            raise SystemExit('published seed derivation not reproduced')
        x_iid = synthetic(g, 'iid', seeds[:200], None)
        x_a = synthetic(g, 'shock', seeds[200:400], g['shock_A'])
        x_mean = mean_stress(g, seeds[400:600])
        x_pub = np.vstack([x_iid, x_a, x_mean])
        obs_pub = observed_stats(x_pub, q, y)
        err = float(np.max(np.abs(obs_pub - saved['observed'])))
        log('published observed reproduction max error', err)
        if err > 1e-8:
            raise SystemExit(f'published observed statistics not reproduced: {err}')
        report['published_observed_max_abs_error'] = err
        un_seeds = {nt: int(np.random.SeedSequence(PUBLISHED_SEED + [tc, k]).generate_state(1)[0]) for nt, k in (('mantel', 71), ('degree', 72))}
    else:
        root = NEW_SEED + [13] if a.context == 'genuine_brain' else NEW_SEED + [15, tc]
        seeds = [int(s.generate_state(1, dtype=np.uint64)[0]) for s in np.random.SeedSequence(root).spawn(400)]
        x_pub = np.vstack([synthetic(g, 'iid', seeds[:200], None), synthetic(g, 'shock', seeds[200:400], g['shock_A'])])
        obs_pub = observed_stats(x_pub, q, y)
        un_seeds = {nt: int(np.random.SeedSequence(NEW_SEED + [14, tc, k]).generate_state(1)[0]) for nt, k in (('mantel', 71), ('degree', 72))}
    b_seeds = [int(s.generate_state(1, dtype=np.uint64)[0]) for s in np.random.SeedSequence(NEW_SEED + [12, tc]).spawn(N_REP)]
    x_b = synthetic(g, 'shock', b_seeds, g['shock_B'])
    obs_b = observed_stats(x_b, q, y)
    xs = np.vstack([x_pub, x_b])
    obs = np.r_[obs_pub, obs_b]
    if a.context in ('brain', 'pbmc'):
        labels = np.array(['iid'] * 200 + ['shock_stratifier_scaled'] * 200 + ['opportunity_mean_stress'] * 200 + ['shock_detection_scaled'] * 200)
    else:
        labels = np.array(['iid'] * 200 + ['shock_stratifier_scaled'] * 200 + ['shock_detection_scaled'] * 200)
    results = {}
    for null in ('unrestricted', 'restricted'):
        nulls = {}
        for nt in ('mantel', 'degree'):
            if null == 'unrestricted':
                gen = null_vectors(g, fixed, od, ind, nt, un_seeds[nt])
            else:
                gen = restricted_vectors(g, fixed, od, ind, nt, NEW_SEED + [11, tc, 0 if nt == 'mantel' else 1])
            log(a.context, null, nt, 'start')
            nulls[nt] = project(xs, fixed, od, g['ii'], nt, gen, a.device)
        if null == 'unrestricted' and a.context in ('brain', 'pbmc'):
            for nt in ('mantel', 'degree'):
                old = np.load(HIST_STRESS / a.context / f'{nt}_nulls.npy')
                e = float(np.max(np.abs(nulls[nt][:, :600] - old)))
                report[f'published_{nt}_null_max_abs_error'] = e
                log('published', nt, 'null reproduction max error', e)
                if e > 1e-8:
                    raise SystemExit(f'published {nt} nulls not reproduced: {e}')
            pm_pub = pvalues(nulls['mantel'][:, :600], obs[:600])
            pd_pub = pvalues(nulls['degree'][:, :600], obs[:600])
            mism = int((pm_pub != saved['pM']).sum() + (pd_pub != saved['pD']).sum())
            report['published_pvalue_mismatches'] = mism
            report['published_pvalues_reproduced_exactly'] = mism == 0
            log('published p-value mismatches', mism)
        pm, pd = pvalues(nulls['mantel'], obs), pvalues(nulls['degree'], obs)
        results[null] = {'pM': pm, 'pD': pd}
        scen = list(dict.fromkeys(labels))
        report[null] = summarize(pm, pd, labels, scen)
        np.savez_compressed(a.out / f'{null}.npz', observed=obs, pM=pm, pD=pd, labels=labels,
                            nullM=nulls['mantel'].astype(np.float32), nullD=nulls['degree'].astype(np.float32))
        (a.out / 'calibration.json').write_text(json.dumps(report, indent=2) + '\n')
    report['status'] = 'complete'
    (a.out / 'calibration.json').write_text(json.dumps(report, indent=2) + '\n')
    log('done', a.context)


if __name__ == '__main__':
    main()
