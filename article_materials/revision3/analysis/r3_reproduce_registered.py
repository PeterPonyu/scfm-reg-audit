#!/usr/bin/env python3
"""Check that the R3 pair-explicit randomization reproduces the registered 999-draw rows.

Brain-selected panel in genuine brain, the four admitted readouts, all non-self pairs,
seed 20261001 per stream: observed effects and both component probabilities must match
the registered `rows.json` of every degree x transformation stream.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3lib as L  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--degree', choices=['none', 'own'], required=True)
    ap.add_argument('--null', choices=['restricted', 'unrestricted'], required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    job = 'brain_selected_panel-brain'
    ctx = L.load_context(job)
    plan = L.load_json(L.VAL / 'corrective_plan.json')
    reg = L.load_json(L.CORR / job / f'statistics-{a.degree}-{a.null}' / 'rows.json')
    names = [r['readout'] for r in reg]
    n = len(ctx['genes'])
    tf = ctx['tf']
    ii, jj = np.repeat(tf, n), np.tile(np.arange(n), len(tf))
    ok = ii != jj
    ii, jj = ii[ok], jj[ok]
    g = ctx['a']['proxy']
    co = ctx['a']['coexpression']
    structural = np.column_stack([ctx['a'][k][jj] for k in ['peak_count', 'gene_length', 'detection', 'gc']])
    controls = np.column_stack([structural, rankdata(co[ii, jj])])
    strata = np.zeros(n, int)
    if a.null == 'restricted':
        ob = np.digitize(ctx['a']['peak_count'], np.unique(np.quantile(ctx['a']['peak_count'], [.2, .4, .6, .8])))
        gb = np.digitize(ctx['a']['gc'], np.unique(np.quantile(ctx['a']['gc'], [1 / 3, 2 / 3])))
        strata = 15 * np.isin(np.arange(n), tf) + 3 * ob + gb
    X = np.stack([ctx['scores'][r][ii, jj] for r in names])
    res = L.components(X, g, tf, ii, jj, controls, strata=strata, permutations=plan['statistics']['permutations'],
                       seed=plan['seed'], designs=[a.degree], tag=f'repro {a.degree}-{a.null}', progress_every=250)
    rows = []
    for k, r in enumerate(reg):
        rows.append({'readout': r['readout'], 'observed_registered': r['observed'], 'observed_r3': float(res['observed'][0, k]),
                     'pM_registered': r['pM'], 'pM_r3': float(res['pM'][0, k]), 'pD_registered': r['pD'], 'pD_r3': float(res['pD'][0, k])})
    ok = all(abs(r['observed_registered'] - r['observed_r3']) < 1e-10 and r['pM_registered'] == r['pM_r3']
             and r['pD_registered'] == r['pD_r3'] for r in rows)
    out = {'job': job, 'degree': a.degree, 'null': a.null, 'permutations': plan['statistics']['permutations'],
           'seed': plan['seed'], 'rows': rows, 'reproduced_exactly': bool(ok), 'basis_impl': L.BASIS_IMPL}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / f'repro-{a.degree}-{a.null}-{L.BASIS_IMPL.replace("+", "_")}.json').write_text(json.dumps(out, indent=1) + '\n')
    print('REPRODUCED' if ok else 'MISMATCH', a.degree, a.null, flush=True)


if __name__ == '__main__':
    main()
