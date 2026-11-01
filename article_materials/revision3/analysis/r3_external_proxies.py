#!/usr/bin/env python3
"""R3-4 inputs: corrective-construction proxies that serve only as external degree sources.

Builds, for the two corrective panels, the proximal motif-accessibility proxy in tissues
other than the tested one, with exactly the corrective rules (strand-aware TSS minus 2 kb
through gene body, 500-bp peak-centred windows, MOODS p < 1e-5, TF OR over motifs and
strands, equal-cell-type log1p mean raw accessibility, zero self-edges):
  brain-selected panel   x GSE206767 fibroblast/induced-neuron pool (pooled; no type labels)
  historical panel       x GSE206767 pool
  historical panel       x PBMC (the registered 1,854 typed PBMC cells)
The existing corrective graphs already supply brain-selected panel x PBMC and both panels x brain.
Only out- and in-degree of these graphs enter any test.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3lib as L  # noqa: E402
import r3paths as P  # noqa: E402
from screg_validation.graphs import read_coordinates, locus_links, build_motif_hits  # noqa: E402
from screg_validation.h5ad import index, read_rows, parse_peak, column  # noqa: E402



def build(panel, tissue, plan, res, selection):
    genes, tfgenes = panel['genes'], panel['tf_genes']
    if tissue == 'fibro':
        apath = P.need('fibroblast_atac')
        with h5py.File(apath) as f:
            cells = index(f, 'obs')
            peaks = index(f, 'var')
            genome = set(column(f, 'var', 'genome'))
        if any('GRCh38' not in v and 'hg38' not in v for v in genome):
            raise SystemExit(f'fibroblast ATAC genome not GRCh38: {genome}')
        asels = list(range(len(cells)))
        labels = ['ALL'] * len(cells)
    else:
        apath = Path(res['pbmc_atac']['path'])
        with h5py.File(apath) as f:
            cells = index(f, 'obs')
            peaks = index(f, 'var')
        look = {c: i for i, c in enumerate(cells)}
        chosen = selection['pbmc_cells']
        asels = [look[r['Barcode']] for r in chosen]
        labels = [r['Cell.Type'] for r in chosen]
    coordinates = read_coordinates(res['gene_coordinates']['path'], plan['panel']['promoter_bp'])
    links = locus_links(peaks, genes, coordinates)
    keep = np.flatnonzero(links.getnnz(1))
    links = links[keep].tocsr()
    windows = []
    w = plan['proxy']['peak_width_bp'] // 2
    for i in keep:
        c, lo, hi = parse_peak(peaks[i])
        mid = (lo + hi) // 2
        windows.append(f'{c}:{max(0, mid - w)}-{mid + w}')
    hits, gc, length = build_motif_hits(windows, tfgenes, res['genome_fasta']['path'], res['meme']['path'],
                                        pvalue=plan['proxy']['motif_pvalue'])
    ax = read_rows(str(apath), 'X', asels, keep)
    if not np.isfinite(ax.data).all() or np.any(ax.data < 0) or not np.allclose(ax.data, np.round(ax.data), atol=1e-6, rtol=0):
        raise SystemExit('ATAC values are not raw integer counts')
    lab = np.array(labels)
    types = sorted(set(labels))
    access = np.mean([np.log1p(np.asarray(ax[lab == t].mean(0)).ravel()) for t in types], axis=0)
    block = (hits.T @ links.multiply(access[:, None])).toarray()
    graph = np.zeros((len(genes), len(genes)))
    graph[[genes.index(g) for g in tfgenes]] = block
    np.fill_diagonal(graph, 0)
    return graph, {'tissue': tissue, 'cells': len(asels), 'types': types, 'linked_peaks': int(len(keep)),
                   'atac_source': str(apath)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--arm', choices=['brain_selected_panel', 'historical_panel_corrected_source'], required=True)
    ap.add_argument('--tissue', choices=['fibro', 'pbmc'], required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    plan = L.load_json(L.VAL / 'corrective_plan.json')
    res = L.load_json(L.VAL / 'corrective_resources.local.json')
    panel = L.load_json(L.CORR / 'panels.json')[a.arm]
    selection = L.load_json(L.CORR / 'selection.json')
    graph, meta = build(panel, a.tissue, plan, res, selection)
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out / f'{a.arm}__{a.tissue}.npz', proxy=graph,
                        out_degree=(graph > 0).sum(1), in_degree=(graph > 0).sum(0))
    L.write_json(a.out / f'{a.arm}__{a.tissue}.json', {**meta, 'arm': a.arm, 'genes': len(panel['genes']),
                 'tfs': len(panel['tf_genes']), 'motif_parameters': plan['proxy'],
                 'scope': 'external degree source only; corrective construction; not a tested outcome'})
    L.log('built', a.arm, a.tissue, meta)


if __name__ == '__main__':
    main()
