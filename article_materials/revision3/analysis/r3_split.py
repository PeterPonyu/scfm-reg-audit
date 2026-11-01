#!/usr/bin/env python3
"""R3-5: disjoint discovery and evaluation samples within GSE174367 (two-fold cross-fitting).

The 18 shared sample identifiers are split into two halves of nine, stratified by
diagnosis x batch with seed 20261031. In fold f, half f is the discovery half: it alone
selects the 1,200-gene panel with the registered rule. The other half is the evaluation
half: its registered RNA nuclei supply the model graphs and co-expression and its
registered ATAC nuclei supply the proxy. Cells are the registered 103 nuclei per sample.

  prepare --fold f   panel from discovery RNA; evaluation RNA/ATAC graphs (CPU)
  infer   --fold f   six readouts on evaluation RNA nuclei (CUDA)
  analyze --fold f --degree d --null n   common support within evaluation half, 9,999 draws
"""
import argparse
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3lib as L  # noqa: E402
from screg_validation.graphs import read_coordinates, read_meme, locus_links, build_motif_hits  # noqa: E402
from screg_validation.h5ad import index, read_rows, parse_peak  # noqa: E402
from screg_validation.models import tokenizable_genes, run_readout  # noqa: E402
from screg_validation.replication import coverage_status  # noqa: E402
from source_correction import read_10x_rows, counts_normalization, subset_genes, valid_coexpression  # noqa: E402

SPLIT_SEED = 20261031
OUT = L.W / 'evidence/split'
JOB = 'split_panel-brain'


def _require_gpu_ignoring_browser(policy):
    """Registered GPU admission, except that a desktop GPU compositor helper process
    (command line containing --type=gpu-process) is not treated as a competing compute job. Temperature, slowdown and recovery checks stay."""
    import os
    import subprocess
    import xml.etree.ElementTree as ET
    import screg_validation.hardware as H
    out = subprocess.run(['nvidia-smi', '-q', '-x'], capture_output=True, text=True, timeout=10)
    root = ET.fromstring(out.stdout)
    ignored = []
    for gpu in root.findall('gpu'):
        procs = gpu.find('processes')
        if procs is None:
            continue
        for p in list(procs):
            name = (p.findtext('process_name') or '')
            if 'type=gpu-process' in name:
                ignored.append({'pid': p.findtext('pid'), 'type': p.findtext('type'), 'name': name.split()[0]})
                procs.remove(p)
    state = H.parse_gpu(ET.tostring(root, encoding='unicode'), policy, own_pid=os.getpid())
    state['ignored_browser_gpu_processes'] = ignored
    if not state['ready']:
        raise RuntimeError('GPU work refused: ' + ', '.join(state['blockers']))
    return state


import screg_validation.models as _models  # noqa: E402
_models.require_gpu = _require_gpu_ignoring_browser


def assignment(records):
    """Stratified halves: within each diagnosis x batch stratum, seeded order, alternate."""
    meta = {}
    for r in records:
        meta.setdefault(r['SampleID'], (r['Diagnosis'], r['Batch']))
    strata = {}
    for s, key in sorted(meta.items()):
        strata.setdefault(key, []).append(s)
    rng = np.random.default_rng(SPLIT_SEED)
    halves = {0: [], 1: []}
    for key in sorted(strata):
        members = list(rng.permutation(sorted(strata[key])))
        for s in members:
            target = 0 if len(halves[0]) < len(halves[1]) else 1 if len(halves[1]) < len(halves[0]) else int(rng.integers(2))
            halves[target].append(str(s))
    if abs(len(halves[0]) - len(halves[1])) > 0:
        raise SystemExit(f'unbalanced halves {halves}')
    return {'half_A': sorted(halves[0]), 'half_B': sorted(halves[1]),
            'strata': {f'{d}|batch{b}': sorted(v) for (d, b), v in strata.items()},
            'diagnosis': {s: meta[s][0] for s in meta}, 'batch': {s: meta[s][1] for s in meta}}


def fold_dir(f):
    return OUT / f'fold{f}'


def prepare(f):
    plan = L.load_json(L.VAL / 'corrective_plan.json')
    res = L.load_json(L.VAL / 'corrective_resources.local.json')
    sel = L.load_json(L.CORR / 'selection.json')
    split = assignment(sel['rna_cells'])
    disc = set(split['half_A'] if f == 1 else split['half_B'])
    evl = set(split['half_B'] if f == 1 else split['half_A'])
    d = fold_dir(f)
    (d / JOB).mkdir(parents=True, exist_ok=True)
    bp = res['brain_rna']['path']
    with h5py.File(bp) as h:
        bc = h['matrix/barcodes'].asstr()[:].tolist()
    look = {c: i for i, c in enumerate(bc)}
    disc_rows = [look[r['Barcode']] for r in sel['rna_cells'] if r['SampleID'] in disc]
    eval_rows = [look[r['Barcode']] for r in sel['rna_cells'] if r['SampleID'] in evl]
    raw_d, allgenes, source = read_10x_rows(bp, disc_rows)
    det = dict(zip(allgenes, np.asarray((raw_d > 0).mean(0)).ravel()))
    coordinates = read_coordinates(res['gene_coordinates']['path'], plan['panel']['promoter_bp'])
    candidates = sorted(set(coordinates) & tokenizable_genes(res))
    with h5py.File(res['brain_atac']['path']) as h:
        ac = index(h, 'obs')
        peaks = index(h, 'var')
    opportunity = locus_links(peaks, candidates, coordinates).getnnz(0)
    detected = [g for g, c in zip(candidates, opportunity) if c > 0 and det.get(g, 0) >= plan['panel']['detection_floor']]
    tf_symbols = {g for m in read_meme(res['meme']['path']).values() for g in m['symbols']}
    tf = sorted(set(detected) & tf_symbols)
    non = sorted(set(detected) - set(tf), key=lambda g: (-det[g], g))
    genes = sorted(tf + non[:plan['panel']['max_genes'] - len(tf)])
    registered = L.load_json(L.CORR / 'panels.json')['brain_selected_panel']
    panel = {'genes': genes, 'tf_genes': tf}
    L.write_json(d / 'panels.json', {'split_panel': panel})
    # Evaluation RNA: full-feature normalization before panel slicing, as registered.
    raw_e, allgenes_e, _ = read_10x_rows(bp, eval_rows)
    cp_e, log_e = counts_normalization(raw_e, 'counts')
    subraw, subcp, sublog = [subset_genes(v, allgenes_e, genes) for v in (raw_e, cp_e, log_e)]
    coexpression, variance = valid_coexpression(sublog)
    de = np.array([float((raw_e[:, allgenes_e.index(g)] > 0).mean()) if g in allgenes_e else 0. for g in genes])
    # Evaluation ATAC nuclei for the proxy.
    alook = {c: i for i, c in enumerate(ac)}
    arecs = [r for r in sel['atac_cells'] if r[plan['brain']['sample_key_atac']] in evl]
    arows = [alook[r['Barcode']] for r in arecs]
    labels = np.array([r['Cell.Type'] for r in arecs])
    links = locus_links(peaks, genes, coordinates)
    keep = np.flatnonzero(links.getnnz(1))
    links = links[keep].tocsr()
    w = plan['proxy']['peak_width_bp'] // 2
    windows = []
    for i in keep:
        c, lo, hi = parse_peak(peaks[i])
        mid = (lo + hi) // 2
        windows.append(f'{c}:{max(0, mid - w)}-{mid + w}')
    hits, gc, length = build_motif_hits(windows, tf, res['genome_fasta']['path'], res['meme']['path'],
                                        pvalue=plan['proxy']['motif_pvalue'])
    ax = read_rows(res['brain_atac']['path'], 'X', arows, keep)
    types = sorted(set(labels))
    access = np.mean([np.log1p(np.asarray(ax[labels == t].mean(0)).ravel()) for t in types], axis=0)
    block = (hits.T @ links.multiply(access[:, None])).toarray()
    graph = np.zeros((len(genes), len(genes)))
    graph[[genes.index(g) for g in tf]] = block
    np.fill_diagonal(graph, 0)
    ng = np.maximum(1, links.getnnz(0))
    gcv = np.asarray(links.T @ gc).ravel() / ng
    sparse.save_npz(d / JOB / 'raw.npz', subraw)
    sparse.save_npz(d / JOB / 'cp10k.npz', subcp)
    np.savez_compressed(d / JOB / 'graphs.npz', proxy=graph, coexpression=coexpression, rna_variance=variance,
                        peak_count=links.getnnz(0), gene_length=np.array([coordinates[g][3] for g in genes]),
                        gc=gcv, detection=de)
    report = {'fold': f, 'discovery_samples': sorted(disc), 'evaluation_samples': sorted(evl), 'split': split,
              'discovery_rna_nuclei': len(disc_rows), 'evaluation_rna_nuclei': len(eval_rows),
              'evaluation_atac_nuclei': len(arows), 'atac_types': types, 'linked_peaks': int(len(keep)),
              'panel_genes': len(genes), 'panel_tfs': len(tf),
              'overlap_with_registered_brain_panel': len(set(genes) & set(registered['genes'])),
              'undefined_coexpression_genes': [g for g, v in zip(genes, variance) if v <= 0],
              'selection_rule': plan['panel']['selection'].replace('genuine brain RNA', 'discovery-half genuine brain RNA'),
              'plan_sha256': (Path(__file__).parent / 'R3_PLAN.sha256').read_text().split()[0]}
    L.write_json(d / 'prepared.json', report)
    L.log('prepared fold', f, {k: report[k] for k in ('discovery_rna_nuclei', 'evaluation_rna_nuclei', 'evaluation_atac_nuclei',
                                                     'panel_genes', 'panel_tfs', 'overlap_with_registered_brain_panel', 'linked_peaks')})


def infer(f):
    plan = L.load_json(L.VAL / 'corrective_plan.json')
    res = L.load_json(L.VAL / 'corrective_resources.local.json')
    d = fold_dir(f)
    panel = L.load_json(d / 'panels.json')['split_panel']
    cp = sparse.load_npz(d / JOB / 'cp10k.npz')
    raw = sparse.load_npz(d / JOB / 'raw.npz')
    for readout in plan['models']['readouts']:
        td = d / JOB / readout
        if (td / 'coverage.json').exists():
            continue
        td.mkdir(exist_ok=True)
        t0 = time.time()
        g, t = run_readout(readout, panel['genes'], cp, raw, res, device='cuda', batch=plan['models']['batch'],
                           seed=plan['seed'], gpu_policy=plan['gpu'])
        status = coverage_status(g, t, panel, plan['models'])
        np.savez_compressed(td / 'graph.npz', graph=g)
        L.write_json(td / 'coverage.json', {**t, **status, 'job': JOB, 'fold': f})
        L.log('inferred fold', f, readout, round(time.time() - t0, 1), 's', status)


def analyze(f, degree, null, permutations):
    d = fold_dir(f)
    ctx = L.load_context(JOB, root=d)
    S, tf_S, ii, jj, info = L.common_support(ctx)
    sub, scores = L.subset(ctx, S)
    X = np.stack([scores[r][ii, jj] for r in L.READOUTS])
    structural = np.column_stack([sub[k][jj] for k in ['peak_count', 'gene_length', 'detection', 'gc']])
    controls = np.column_stack([structural, rankdata(sub['coexpression'][ii, jj])])
    strata = L.strata_for(sub, tf_S, null == 'restricted')
    t0 = time.time()
    r = L.components(X, sub['proxy'], tf_S, ii, jj, controls, strata=strata, permutations=permutations,
                     seed=L.SEED, designs=[degree], tag=f'split fold{f} {degree}-{null}', progress_every=1000)
    rows = [{'readout': n, 'observed': float(r['observed'][0, k]), 'pM': float(r['pM'][0, k]),
             'pD': float(r['pD'][0, k]), 'pIUT': float(r['pIUT'][0, k])} for k, n in enumerate(L.READOUTS)]
    (L.W / 'evidence/split_results').mkdir(parents=True, exist_ok=True)
    L.write_json(L.W / 'evidence/split_results' / f'split__fold{f}__{degree}-{null}.json',
                 {'fold': f, 'degree_mode': degree, 'null_mode': null, 'permutations': permutations, 'seed': L.SEED,
                  'support': info, 'rows': rows, 'seconds': round(time.time() - t0, 1), 'basis_impl': L.BASIS_IMPL})
    L.log('analyzed fold', f, degree, null)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('task', choices=['prepare', 'infer', 'analyze'])
    ap.add_argument('--fold', type=int, choices=[1, 2], required=True)
    ap.add_argument('--degree', choices=['none', 'own'])
    ap.add_argument('--null', choices=['restricted', 'unrestricted'])
    ap.add_argument('--permutations', type=int, default=9999)
    a = ap.parse_args()
    if a.task == 'prepare':
        prepare(a.fold)
    elif a.task == 'infer':
        infer(a.fold)
    else:
        analyze(a.fold, a.degree, a.null, a.permutations)


if __name__ == '__main__':
    main()
