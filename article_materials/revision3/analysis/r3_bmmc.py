#!/usr/bin/env python3
"""R3-6: five-donor bone-marrow analysis on the common finite gene set.

Gene set: panel genes with positive RNA variance and finite gene-level scores in all five
evaluation donors for the five gene-level readouts (Geneformer attention stays unavailable:
its pair coverage is 73-93% per donor). Degree modes reference/own/none, restricted and
unrestricted components, 9,999 draws, seed 20260930 per stream, as registered. Reference
degree is the mean of the two degree-reference donors' proxy degrees within the same gene set.

  run   --donor D --degree d --null n
  check --draws K   reproduce the first K registered donor-16710 draws (full panel)
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3lib as L  # noqa: E402

REP = L.VAL / 'runs/a5c052c44a1d1ce7ba63/replication'
FIVE = ['geneformer_embed', 'random_init_geneformer', 'scgpt_encoder', 'scfoundation_encoder', 'uce_encoder']
OUT = L.W / 'evidence/bmmc'


def setup():
    plan = L.load_json(L.VAL / 'plan.json')
    panel = L.load_json(REP / 'panel.json')
    prepared = L.load_json(REP / 'prepared.json')
    conf = dict(np.load(REP / 'confounds.npz', allow_pickle=False))
    return plan, panel, prepared, conf


def donor_arrays(prepared, donor):
    stem = prepared['donor_stems'][donor]
    a = dict(np.load(REP / f'{stem}-graphs.npz', allow_pickle=False))
    cp = sparse.load_npz(REP / f'{stem}-cp10k.npz')
    lg = cp.copy()
    lg.data = np.log1p(lg.data)
    var = np.asarray(lg.multiply(lg).mean(0)).ravel() - np.asarray(lg.mean(0)).ravel() ** 2
    scores = {r: np.load(REP / f'{stem}-{r}' / 'model_graph.npz', allow_pickle=False)['graph'] for r in FIVE}
    return a, var, scores


def common_genes(plan, panel, prepared):
    n = len(panel['genes'])
    keep = np.ones(n, bool)
    detail = {}
    for d in panel['split']['evaluation']:
        a, var, scores = donor_arrays(prepared, d)
        k = var > 1e-12
        for r in FIVE:
            k &= ~L.missing_genes(scores[r])
        detail[d] = [panel['genes'][i] for i in np.flatnonzero(~k)]
        keep &= k
    return np.flatnonzero(keep), detail


def restricted_strata(conf, S, tf_S):
    m = len(S)
    ob = np.digitize(conf['opportunity_degree'][S], np.unique(np.quantile(conf['opportunity_degree'][S], [.2, .4, .6, .8])))
    gb = np.digitize(conf['gc'][S], np.unique(np.quantile(conf['gc'][S], [1 / 3, 2 / 3])))
    return ob * 3 + gb + 15 * np.isin(np.arange(m), tf_S)


def pairs_for(tf_idx, m):
    ii, jj = np.repeat(tf_idx, m), np.tile(np.arange(m), len(tf_idx))
    ok = ii != jj
    return ii[ok], jj[ok]


def run(donor, degree, null, permutations):
    plan, panel, prepared, conf = setup()
    S, detail = common_genes(plan, panel, prepared)
    pos = -np.ones(len(panel['genes']), int)
    pos[S] = np.arange(len(S))
    tf_full = [panel['genes'].index(t) for t in panel['tf_genes']]
    tf_S = np.array([pos[t] for t in tf_full if pos[t] >= 0])
    ii, jj = pairs_for(tf_S, len(S))
    a, var, scores = donor_arrays(prepared, donor)
    g = a['proxy'][np.ix_(S, S)]
    co = a['coexpression'][np.ix_(S, S)]
    common = np.column_stack([conf[k][S][jj] for k in ('peak_count', 'gene_length', 'detection', 'gc')])
    controls = np.column_stack([common, rankdata(co[ii, jj])])
    refs = []
    for r in panel['split']['degree_reference']:
        gr = np.load(REP / f"{prepared['donor_stems'][r]}-graphs.npz", allow_pickle=False)['proxy'][np.ix_(S, S)]
        refs.append(((gr > 0).sum(1), (gr > 0).sum(0)))
    od_ref = np.mean([x[0] for x in refs], axis=0)
    ind_ref = np.mean([x[1] for x in refs], axis=0)
    strata = restricted_strata(conf, S, tf_S) if null == 'restricted' else np.zeros(len(S), int)
    X = np.stack([scores[r][np.ix_(S, S)][ii, jj] for r in FIVE])
    backend = 'cuda' if (Path(__file__).parent / 'USE_CUDA').exists() else 'cpu'
    t0 = time.time()
    res = L.components(X, g, tf_S, ii, jj, controls, strata=strata, permutations=permutations, seed=plan['seed'],
                       designs=[degree], external={'reference': (od_ref, ind_ref)},
                       tag=f'bmmc {donor} {degree}-{null}', progress_every=1000, backend=backend)
    base_c = np.column_stack([common, od_ref[ii], ind_ref[jj]])
    base_rho = L.partial_spearman(co[ii, jj], g[ii, jj], base_c)
    rows = []
    for k, r in enumerate(FIVE):
        fm = L.partial_spearman(X[k], g[ii, jj], base_c)
        rows.append({'donor': donor, 'readout': r, 'degree_mode': degree, 'null_mode': null,
                     'observed': float(res['observed'][0, k]), 'pM': float(res['pM'][0, k]), 'pD': float(res['pD'][0, k]),
                     'pIUT': float(res['pIUT'][0, k]), 'comparable_marginal_fm': fm, 'comparable_marginal_baseline': base_rho,
                     'comparable_rho_difference': fm - base_rho,
                     'concordance': float(spearmanr(X[k], co[ii, jj]).statistic)})
    OUT.mkdir(parents=True, exist_ok=True)
    L.write_json(OUT / f'bmmc__{donor}__{degree}-{null}.json', {
        'donor': donor, 'degree_mode': degree, 'null_mode': null, 'permutations': permutations, 'seed': plan['seed'],
        'genes_common': int(len(S)), 'tfs_common': int(len(tf_S)), 'pairs': int(len(ii)),
        'genes_panel': len(panel['genes']), 'tfs_panel': len(panel['tf_genes']), 'excluded_by_donor': detail,
        'rows': rows, 'seconds': round(time.time() - t0, 1), 'basis_impl': L.BASIS_IMPL, 'backend': backend})
    L.log('done bmmc', donor, degree, null)


def check(draws):
    """Registered donor 16710 used the full panel; compare our first draws with its checkpoints."""
    plan, panel, prepared, conf = setup()
    donor = '16710'
    n = len(panel['genes'])
    tf = np.array([panel['genes'].index(t) for t in panel['tf_genes']])
    ii, jj = pairs_for(tf, n)
    a, var, scores = donor_arrays(prepared, donor)
    common = np.column_stack([conf[k][jj] for k in ('peak_count', 'gene_length', 'detection', 'gc')])
    controls = np.column_stack([common, rankdata(a['coexpression'][ii, jj])])
    reference = [np.load(REP / f"{prepared['donor_stems'][r]}-graphs.npz", allow_pickle=False)['proxy']
                 for r in panel['split']['degree_reference']]
    od = np.mean([(g > 0).sum(1) for g in reference], axis=0)
    ind = np.mean([(g > 0).sum(0) for g in reference], axis=0)
    report = {}
    for degree, null in (('reference', 'restricted'), ('own', 'unrestricted')):
        stem = prepared['donor_stems'][donor]
        state = L.load_json(REP / f'{stem}-statistics-{degree}-{null}' / 'draw-progress.json')
        rows = L.load_json(REP / f'{stem}-statistics-{degree}-{null}' / 'rows.json')
        names = [r['readout'] for r in rows]
        conf_full = {'peak_count': conf['peak_count'], 'gc': conf['gc'], 'opportunity_degree': conf['opportunity_degree']}
        strata = restricted_strata(conf_full, np.arange(n), tf) if null == 'restricted' else np.zeros(n, int)
        X = np.stack([np.load(REP / f'{stem}-{r}' / 'model_graph.npz', allow_pickle=False)['graph'][ii, jj] for r in names])
        res = L.components(X, a['proxy'], tf, ii, jj, controls, strata=strata, permutations=draws, seed=plan['seed'],
                           designs=[degree], external={'reference': (od, ind)}, tag=f'check {degree}-{null}', progress_every=0)
        savedM = np.asarray(state['nulls']['M'], float).reshape(state['next_draw'], len(names))[:draws]
        savedD = np.asarray(state['nulls']['D'], float).reshape(state['next_draw'], len(names))[:draws]
        eM = float(np.max(np.abs(res['nulls']['M'][:, 0, :] - savedM)))
        eD = float(np.max(np.abs(res['nulls']['D'][:, 0, :] - savedD)))
        eo = float(np.max(np.abs(res['observed'][0] - np.array([r['observed'] for r in rows]))))
        report[f'{degree}-{null}'] = {'observed_max_abs_error': eo, 'nullM_max_abs_error': eM, 'nullD_max_abs_error': eD,
                                      'draws_compared': draws}
        L.log('check', degree, null, report[f'{degree}-{null}'])
    OUT.mkdir(parents=True, exist_ok=True)
    L.write_json(OUT / 'reproduction_16710.json', report)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('task', choices=['run', 'check', 'genes'])
    ap.add_argument('--donor')
    ap.add_argument('--degree', choices=['reference', 'own', 'none'])
    ap.add_argument('--null', choices=['restricted', 'unrestricted'])
    ap.add_argument('--permutations', type=int, default=9999)
    ap.add_argument('--draws', type=int, default=256)
    a = ap.parse_args()
    if a.task == 'run':
        run(a.donor, a.degree, a.null, a.permutations)
    elif a.task == 'check':
        check(a.draws)
    else:
        plan, panel, prepared, conf = setup()
        S, detail = common_genes(plan, panel, prepared)
        tf_S = [t for t in panel['tf_genes'] if panel['genes'].index(t) in set(S)]
        print('common genes', len(S), 'TFs', len(tf_S), {d: len(v) for d, v in detail.items()})
        print('excluded TFs', sorted(set(panel['tf_genes']) - set(tf_S)))


if __name__ == '__main__':
    main()
