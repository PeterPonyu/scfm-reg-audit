#!/usr/bin/env python3
"""R3-2 / R3-3 streams: common-support family, attention-free sensitivity, symmetric baseline.

variant six        : six readouts on pairs E (identical edge sets), controls = structural + ranked co-expression
variant gene_level : five gene-level readouts on all pairs within S, same controls
variant symmetric  : six readouts plus co-expression as a seventh score, controls = structural only
One stream = one job x degree mode x transformation class; 9,999 draws, seed 20261001.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3lib as L  # noqa: E402

STRUCT = ['peak_count', 'gene_length', 'detection', 'gc']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--job', choices=L.JOBS, required=True)
    ap.add_argument('--degree', choices=['none', 'own'], required=True)
    ap.add_argument('--null', choices=['restricted', 'unrestricted'], required=True)
    ap.add_argument('--variant', choices=['six', 'gene_level', 'symmetric'], required=True)
    ap.add_argument('--permutations', type=int, default=9999)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    t0 = time.time()
    ctx = L.load_context(a.job)
    pair_readouts = L.GENE_LEVEL if a.variant == 'gene_level' else L.READOUTS
    S, tf_S, ii, jj, info = L.common_support(ctx, readouts=L.READOUTS, pair_readouts=pair_readouts)
    sub, scores = L.subset(ctx, S)
    names = list(L.GENE_LEVEL if a.variant == 'gene_level' else L.READOUTS)
    X = [scores[r][ii, jj] for r in names]
    structural = np.column_stack([sub[k][jj] for k in STRUCT])
    co = sub['coexpression'][ii, jj]
    if a.variant == 'symmetric':
        controls = structural
        X.append(co)
        names.append('coexpression')
    else:
        controls = np.column_stack([structural, rankdata(co)])
    strata = L.strata_for(sub, tf_S, a.null == 'restricted')
    backend = 'cuda' if (Path(__file__).parent / 'USE_CUDA').exists() else 'cpu'
    res = L.components(np.stack(X), sub['proxy'], tf_S, ii, jj, controls, strata=strata,
                       permutations=a.permutations, seed=L.SEED, designs=[a.degree],
                       tag=f'{a.job} {a.variant} {a.degree}-{a.null}', progress_every=1000, backend=backend)
    rows = []
    for k, r in enumerate(names):
        rows.append({'readout': r, 'observed': float(res['observed'][0, k]), 'pM': float(res['pM'][0, k]),
                     'pD': float(res['pD'][0, k]), 'pIUT': float(res['pIUT'][0, k])})
    a.out.mkdir(parents=True, exist_ok=True)
    stem = f'{a.variant}__{a.job}__{a.degree}-{a.null}'
    L.write_json(a.out / f'{stem}.json', {
        'job': a.job, 'context': L.LABEL[a.job], 'variant': a.variant, 'degree_mode': a.degree, 'null_mode': a.null,
        'permutations': a.permutations, 'seed': L.SEED, 'support': info, 'stratum_sizes': res['stratum_sizes'],
        'controls': (STRUCT if a.variant == 'symmetric' else STRUCT + ['ranked co-expression']),
        'rows': rows, 'seconds': round(time.time() - t0, 1), 'basis_impl': L.BASIS_IMPL, 'backend': backend,
        'plan_sha256': (Path(__file__).parent / 'R3_PLAN.sha256').read_text().split()[0]})
    np.savez_compressed(a.out / f'{stem}__nulls.npz', M=res['nulls']['M'][:, 0, :].astype(np.float32),
                        D=res['nulls']['D'][:, 0, :].astype(np.float32), names=np.array(names))
    L.log('done', stem, round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
