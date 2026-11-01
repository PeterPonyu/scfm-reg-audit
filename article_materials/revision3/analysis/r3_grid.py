#!/usr/bin/env python3
"""R3-4: the ten historical degree designs on the three source-matched corrective contexts.

Pairs E (identical edge sets for all six readouts), controls = structural + ranked
co-expression, unrestricted components as in the historical grid, 9,999 draws, seed
20261001. Designs can be split across processes: the random stream depends only on the
seed and strata, so every process sees identical draws.
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
EXTERNAL = {
    'brain_selected_panel-brain': {'external_fibro': ('ext', 'brain_selected_panel__fibro'),
                                   'external_other': ('corr', 'brain_selected_panel-pbmc')},
    'brain_selected_panel-pbmc': {'external_fibro': ('ext', 'brain_selected_panel__fibro'),
                                  'external_other': ('corr', 'brain_selected_panel-brain')},
    'historical_panel_corrected_source-brain': {'external_fibro': ('ext', 'historical_panel_corrected_source__fibro'),
                                                'external_other': ('ext', 'historical_panel_corrected_source__pbmc')},
}


def external_degrees(job, S):
    out = {}
    for name, (kind, stem) in EXTERNAL[job].items():
        if kind == 'ext':
            g = np.load(L.W / 'evidence/external_proxies' / f'{stem}.npz', allow_pickle=False)['proxy']
        else:
            g = np.load(L.CORR / stem / 'graphs.npz', allow_pickle=False)['proxy']
        od, ind = (g > 0).sum(1).astype(float), (g > 0).sum(0).astype(float)
        out[name] = (od[S], ind[S])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--job', choices=L.JOBS, required=True)
    ap.add_argument('--designs', nargs='+', choices=L.DESIGNS, required=True)
    ap.add_argument('--permutations', type=int, default=9999)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
    a = ap.parse_args()
    t0 = time.time()
    ctx = L.load_context(a.job)
    S, tf_S, ii, jj, info = L.common_support(ctx)
    sub, scores = L.subset(ctx, S)
    X = np.stack([scores[r][ii, jj] for r in L.READOUTS])
    controls = np.column_stack([np.column_stack([sub[k][jj] for k in STRUCT]), rankdata(sub['coexpression'][ii, jj])])
    ext = external_degrees(a.job, S)
    strata = np.zeros(len(S), int)
    res = L.components(X, sub['proxy'], tf_S, ii, jj, controls, strata=strata, permutations=a.permutations,
                       seed=L.SEED, designs=a.designs, external=ext, tag=f'grid {a.job} {"+".join(a.designs)}',
                       progress_every=1000, backend=a.backend)
    rows = []
    for di, d in enumerate(a.designs):
        for k, r in enumerate(L.READOUTS):
            rows.append({'design': d, 'readout': r, 'observed': float(res['observed'][di, k]),
                         'pM': float(res['pM'][di, k]), 'pD': float(res['pD'][di, k]), 'pIUT': float(res['pIUT'][di, k]),
                         'nullM_sd': float(res['nulls']['M'][:, di, k].std()), 'nullD_sd': float(res['nulls']['D'][:, di, k].std())})
    a.out.mkdir(parents=True, exist_ok=True)
    stem = f'grid__{a.job}__{"+".join(a.designs)}'
    L.write_json(a.out / f'{stem}.json', {'job': a.job, 'context': L.LABEL[a.job], 'designs': a.designs,
                 'null_mode': 'unrestricted', 'permutations': a.permutations, 'seed': L.SEED, 'support': info,
                 'external_degree_sources': {k: v[1] for k, v in EXTERNAL[a.job].items()}, 'rows': rows, 'basis_impl': L.BASIS_IMPL, 'backend': a.backend,
                 'seconds': round(time.time() - t0, 1),
                 'plan_sha256': (Path(__file__).parent / 'R3_PLAN.sha256').read_text().split()[0]})
    L.log('done', stem, round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
