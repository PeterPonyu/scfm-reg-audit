#!/usr/bin/env python3
"""Aggregate every R3 stream into the families, corrections and summaries used in the paper.

Writes evidence/R3_RESULTS.json. Families are exactly those fixed in R3_PLAN.md; every
slot is measured. BY and Holm use the registered `adjust` implementation.
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3lib as L  # noqa: E402
from screg_validation.statistics import adjust, donor_interval, exact_paired_signflip  # noqa: E402

E = L.W / 'evidence'
STRUCT = ['peak_count', 'gene_length', 'detection', 'gc']


def fam(rows, key='pIUT'):
    p = [r[key] for r in rows]
    by, holm = adjust(p), adjust(p, 'Holm')
    for r, q, h in zip(rows, by, holm):
        r['q_BY'], r['p_Holm'] = float(q), float(h)
    sel = [r for r in rows if r['q_BY'] < .05 or r['p_Holm'] < .05]
    return {'slots': len(rows), 'min_raw_pIUT': float(min(p)), 'min_BY': float(min(by)), 'min_Holm': float(min(holm)),
            'BY_selected': int(sum(q < .05 for q in by)), 'Holm_selected': int(sum(h < .05 for h in holm)),
            'selected': [{k: r[k] for k in r if k in ('context', 'job', 'readout', 'degree_mode', 'null_mode', 'design',
                                                      'fold', 'donor', 'observed', 'pIUT', 'q_BY', 'p_Holm')} for r in sel],
            'holm_floor': float(len(rows) / 10000)}


def load_streams(pattern):
    out = []
    for f in sorted(glob.glob(str(pattern))):
        if f.endswith('__nulls.npz'):
            continue
        out.append(json.loads(Path(f).read_text()))
    return out


def common_support():
    res = {}
    streams = load_streams(E / 'common_support' / '*.json')
    for variant, expect in (('six', 12), ('gene_level', 12), ('symmetric', 12)):
        ss = [s for s in streams if s['variant'] == variant]
        rows = []
        for s in ss:
            for r in s['rows']:
                if variant == 'symmetric' and r['readout'] == 'coexpression':
                    continue
                rows.append({'variant': variant, 'job': s['job'], 'context': s['context'], 'degree_mode': s['degree_mode'],
                             'null_mode': s['null_mode'], **r})
        entry = {'streams_found': len(ss), 'streams_expected': expect, 'rows': rows,
                 'support': {s['job']: s['support'] for s in ss}}
        if len(ss) == expect:
            entry['family'] = fam(rows)
            entry['family_restricted_only'] = fam([dict(r) for r in rows if r['null_mode'] == 'restricted'])
            for r in rows:  # restore 72-slot adjusted values on the main rows
                pass
        if variant == 'symmetric':
            entry['coexpression_rows'] = [{'job': s['job'], 'context': s['context'], 'degree_mode': s['degree_mode'],
                                           'null_mode': s['null_mode'], **r} for s in ss for r in s['rows']
                                          if r['readout'] == 'coexpression']
        res[variant] = entry
    if res['six'].get('family'):
        fam(res['six']['rows'])  # 72-slot values are the reported ones
    if res['gene_level'].get('family'):
        fam(res['gene_level']['rows'])
    if res['symmetric'].get('family'):
        fam(res['symmetric']['rows'])
    return res


def tf_bootstrap(job, degree, resamples=2000, seed=20261031):
    """Symmetric contrast rho(FM|C) - rho(coexpression|C), C = structural (+ observed own degree).
    TF-cluster bootstrap of residualized rank correlations (residualization fixed at the full fit)."""
    ctx = L.load_context(job)
    S, tf_S, ii, jj, info = L.common_support(ctx)
    sub, scores = L.subset(ctx, S)
    g = sub['proxy']
    C = np.column_stack([np.ones(len(ii))] + [sub[k][jj] for k in STRUCT])
    if degree == 'own':
        C = np.column_stack([C, (g > 0).sum(1)[ii], (g > 0).sum(0)[jj]])
    q = L.basis(C)

    def res_rank(v):
        r = rankdata(v)
        r = r - r.mean()
        return r - q @ (q.T @ r)
    y = res_rank(g[ii, jj])
    c = res_rank(sub['coexpression'][ii, jj])
    rows = np.unique(ii)
    inv = np.searchsorted(rows, ii)
    T = len(rows)

    def sums(a, b):
        return np.bincount(inv, weights=a * b, minlength=T)
    n_t = np.bincount(inv, minlength=T).astype(float)
    Sy, Sc = np.bincount(inv, weights=y, minlength=T), np.bincount(inv, weights=c, minlength=T)
    Syy, Scc, Scy = sums(y, y), sums(c, c), sums(c, y)
    rng = np.random.default_rng(seed)
    W = np.stack([np.bincount(rng.integers(0, T, size=T), minlength=T) for _ in range(resamples)]).astype(float)

    def corr(Sab, Sa, Sb, Saa, Sbb):
        n = W @ n_t
        cov = W @ Sab - (W @ Sa) * (W @ Sb) / n
        va = W @ Saa - (W @ Sa) ** 2 / n
        vb = W @ Sbb - (W @ Sb) ** 2 / n
        return cov / np.sqrt(va * vb)
    rho_c_obs = float(c @ y / np.sqrt((c @ c) * (y @ y)))
    boot_c = corr(Scy, Sc, Sy, Scc, Syy)
    out = {'job': job, 'context': L.LABEL[job], 'degree_mode': degree, 'coexpression_rho': rho_c_obs,
           'coexpression_ci95': [float(v) for v in np.quantile(boot_c, [.025, .975])], 'readouts': []}
    for r in L.READOUTS:
        x = res_rank(scores[r][ii, jj])
        Sx, Sxx, Sxy = np.bincount(inv, weights=x, minlength=T), sums(x, x), sums(x, y)
        rho = float(x @ y / np.sqrt((x @ x) * (y @ y)))
        boot_x = corr(Sxy, Sx, Sy, Sxx, Syy)
        d = boot_x - boot_c
        out['readouts'].append({'readout': r, 'fm_rho': rho, 'fm_ci95': [float(v) for v in np.quantile(boot_x, [.025, .975])],
                                'difference': rho - rho_c_obs,
                                'difference_ci95': [float(v) for v in np.quantile(d, [.025, .975])]})
    out['resamples'] = resamples
    out['tfs_resampled'] = int(T)
    return out


def grid():
    streams = load_streams(E / 'grid' / '*.json')
    rows = []
    for s in streams:
        for r in s['rows']:
            rows.append({'job': s['job'], 'context': s['context'], **r})
    designs_found = {(r['job'], r['design']) for r in rows}
    entry = {'streams_found': len(streams), 'rows': rows, 'designs_found': len(designs_found)}
    if len(designs_found) == 30:
        per = {}
        for job in L.JOBS:
            for d in L.DESIGNS:
                sub = [dict(r) for r in rows if r['job'] == job and r['design'] == d]
                f = fam(sub)
                per[f'{job}|{d}'] = {'BY_selected': f['BY_selected'], 'Holm_selected': f['Holm_selected'],
                                     'min_BY': f['min_BY'], 'min_Holm': f['min_Holm']}
        entry['per_context_design'] = per
        entry['selection_by_design'] = {d: {'BY': sum(per[f'{j}|{d}']['BY_selected'] for j in L.JOBS),
                                            'Holm': sum(per[f'{j}|{d}']['Holm_selected'] for j in L.JOBS)} for d in L.DESIGNS}
        entry['pooled_180'] = fam(rows)
    return entry


def split():
    streams = load_streams(E / 'split_results' / '*.json')
    rows = [{'fold': s['fold'], 'degree_mode': s['degree_mode'], 'null_mode': s['null_mode'], **r}
            for s in streams for r in s['rows']]
    prep = {f: json.loads((E / 'split' / f'fold{f}' / 'prepared.json').read_text()) for f in (1, 2)
            if (E / 'split' / f'fold{f}' / 'prepared.json').exists()}
    entry = {'streams_found': len(streams), 'rows': rows, 'support': {s['fold']: s['support'] for s in streams},
             'prepared': {f: {k: v for k, v in p.items() if k != 'split'} for f, p in prep.items()},
             'assignment': prep.get(1, {}).get('split')}
    if len(streams) == 8:
        entry['family_48'] = fam(rows)
        entry['family_restricted_24'] = fam([dict(r) for r in rows if r['null_mode'] == 'restricted'])
        fam(rows)
    return entry


def bmmc():
    streams = load_streams(E / 'bmmc' / 'bmmc__*.json')
    rows = [r for s in streams for r in s['rows']]
    entry = {'streams_found': len(streams), 'rows': rows}
    if streams:
        entry['common'] = {k: streams[0][k] for k in ('genes_common', 'tfs_common', 'pairs', 'genes_panel', 'tfs_panel', 'excluded_by_donor')}
    if len(streams) == 30:
        entry['family_150'] = fam(rows)
        entry['family_restricted_75'] = fam([dict(r) for r in rows if r['null_mode'] == 'restricted'])
        fam(rows)
        donors = sorted({r['donor'] for r in rows})
        summary = []
        for readout in sorted({r['readout'] for r in rows}):
            base = [r for r in rows if r['readout'] == readout and r['degree_mode'] == 'reference' and r['null_mode'] == 'restricted']
            base = sorted(base, key=lambda r: r['donor'])
            vals = [r['observed'] for r in base]
            diffs = [r['comparable_rho_difference'] for r in base]
            summary.append({'readout': readout, 'donors': [r['donor'] for r in base], 'observed': vals,
                            'donor_interval': donor_interval(vals, [r['donor'] for r in base], expected_donors=donors, resamples=2000, seed=20260930),
                            'differences': diffs,
                            'paired_baseline_interval': donor_interval(diffs, [r['donor'] for r in base], expected_donors=donors, resamples=2000, seed=20260930),
                            'paired_baseline_signflip_p': exact_paired_signflip(diffs)})
        holm = adjust([s['paired_baseline_signflip_p'] for s in summary], 'Holm')
        for s, h in zip(summary, holm):
            s['paired_baseline_signflip_Holm'] = float(h)
        entry['donor_summary'] = summary
        entry['minimum_two_sided_signflip_p'] = 2 / 2 ** 5
    return entry


def main():
    out = {'plan_sha256': (L.W / 'analysis/R3_PLAN.sha256').read_text().split()[0]}
    out['calibration'] = {c: json.loads((E / 'calibration' / c / 'calibration.json').read_text())
                          for c in ('brain', 'pbmc', 'genuine_brain', 'new_pbmc', 'old_brain', 'bmmc_13272', 'bmmc_28483')
                          if (E / 'calibration' / c / 'calibration.json').exists()}
    out['reproduction'] = {'registered_brain_streams': [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(E / 'reproduction/*.json')))],
                           'bmmc_16710': json.loads((E / 'bmmc/reproduction_16710.json').read_text()) if (E / 'bmmc/reproduction_16710.json').exists() else None}
    out['common_support'] = common_support()
    if '--bootstrap' in sys.argv:
        out['symmetric_bootstrap'] = [tf_bootstrap(j, d) for j in L.JOBS for d in ('none', 'own')]
    elif (E / 'R3_RESULTS.json').exists():
        prev = json.loads((E / 'R3_RESULTS.json').read_text())
        if 'symmetric_bootstrap' in prev:
            out['symmetric_bootstrap'] = prev['symmetric_bootstrap']
    out['grid'] = grid()
    out['split'] = split()
    out['bmmc'] = bmmc()
    (E / 'R3_RESULTS.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    cs = out['common_support']
    print('common support streams', {k: v['streams_found'] for k, v in cs.items()}, '| grid', out['grid']['designs_found'],
          '| split', out['split']['streams_found'], '| bmmc', out['bmmc']['streams_found'])


if __name__ == '__main__':
    main()
