#!/usr/bin/env python3
"""Official sample counts identity audit with exact feature/barcode alignment."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
from scipy import sparse
from screg_validation.common import ContractError, digest, load_json, sha256, write_json
from screg_validation.h5ad import index, column, read_rows


def unique_symbols(symbols):
    taken = set(symbols)
    seen = set()
    counts = {}
    result = []
    for value in symbols:
        if value not in seen:
            result.append(value)
            seen.add(value)
            continue
        i = counts.get(value, 0) + 1
        while value + '-' + str(i) in taken:
            i += 1
        renamed = value + '-' + str(i)
        counts[value] = i
        taken.add(renamed)
        result.append(renamed)
    return result


def counts_fingerprint(matrix_path, barcode_path, features_path, candidate, *, preceding_barcodes, max_cells=128, seed=20260930):
    for path in (matrix_path, barcode_path, features_path, preceding_barcodes):
        record = load_json(Path(str(path) + '.download.json'))
        if record['sha256'] != sha256(path) or record['bytes'] != Path(path).stat().st_size:
            raise ContractError('official source download identity failed')
    with gzip.open(barcode_path, 'rt') as f:
        barcodes = [l.strip() for l in f]
    with gzip.open(features_path, 'rt') as f:
        features = [l.strip().split('\t') for l in f]
    with h5py.File(candidate, 'r') as f:
        candidate_cells = index(f, 'obs')
        candidate_genes = index(f, 'var')
        batches = column(f, 'obs', 'batch') if 'batch' in f['obs'] else None
    official_genes = unique_symbols([row[1] for row in features])
    if official_genes != candidate_genes or len(barcodes) != len(set(barcodes)):
        raise ContractError('official feature order or barcode uniqueness differs')
    with gzip.open(preceding_barcodes, 'rt') as f:
        prior_barcodes = [l.strip() for l in f]
    # Reconstruct the entire observed Normal1 then AD1 concatenation. This is
    # exact file-backed naming replay, not generic stripping of barcode suffixes.
    combined = unique_symbols(prior_barcodes + barcodes)
    if combined != candidate_cells:
        raise ContractError('complete official concatenation does not replay candidate cell order')
    offset = len(prior_barcodes)
    mapping = {barcodes[i]: (offset + i, candidate_cells[offset + i]) for i in range(len(barcodes))}
    common = sorted(barcodes,
                    key=lambda c: hashlib.sha256(f'{seed}:GSM6593315:{c}'.encode()).hexdigest())[:max_cells]
    official_lookup = {c:i for i,c in enumerate(barcodes)}
    candidate_lookup = {c:mapping[c][0] for c in common}
    selected_columns = {official_lookup[c]:i for i,c in enumerate(common)}
    rows, cols, data = [], [], []
    dimensions = None
    with gzip.open(matrix_path, 'rt') as f:
        header = f.readline().strip()
        if header != '%%MatrixMarket matrix coordinate integer general':
            raise ContractError('explicit integer MatrixMarket source required')
        for line in f:
            if line.startswith('%'):
                continue
            dimensions = tuple(map(int, line.split()))
            break
        if dimensions is None or dimensions[:2] != (len(features),len(barcodes)):
            raise ContractError('official feature/barcode matrix dimensions differ')
        entries = 0
        for line in f:
            if not line.strip():continue
            rr, cc, vv = map(int, line.split());entries += 1
            if not 1 <= rr <= len(features) or not 1 <= cc <= len(barcodes) or vv < 0:
                raise ContractError('invalid raw count coordinate')
            if cc-1 in selected_columns:
                rows.append(selected_columns[cc-1]);cols.append(rr-1);data.append(vv)
        if entries != dimensions[2]:
            raise ContractError('MatrixMarket declared nonzero count differs')
    official = sparse.csr_matrix((data,(rows,cols)),shape=(len(common),len(features)))
    local = read_rows(candidate,'X',[candidate_lookup[c] for c in common])
    diff=(official-local).tocsr();diff.eliminate_zeros()
    matches = []
    for i,c in enumerate(common):
        dr=diff[i];matches.append({'official_sample':'GSM6593315','official_barcode':c,'local_cell':mapping[c][1],'exact_counts_match':dr.nnz==0,'max_absolute_difference':float(abs(dr.data).max()) if dr.nnz else 0.,'local_batch':batches[candidate_lookup[c]] if batches else None})
    return {'status':'official_counts_match' if diff.nnz==0 else 'official_counts_conflict','official_matrix_sha256':sha256(matrix_path),'official_barcodes_sha256':sha256(barcode_path),'official_features_sha256':sha256(features_path),'candidate_sha256':sha256(candidate),'dimensions':list(dimensions),'matched_identity_cells':len(barcodes),'naming_replay':'official Normal1 then AD1 concatenation with exact make_index_unique replay across all candidate rows','preceding_barcodes_sha256':sha256(preceding_barcodes),'cells_compared':len(common),'all_feature_count':len(features),'exact_numeric_cells':sum(r['exact_counts_match'] for r in matches),'max_absolute_difference':float(abs(diff.data).max()) if diff.nnz else 0.,'rows':matches,'scope':'Full raw feature vector identity for deterministic bounded cells in GSM6593315, after exact composite sample-barcode naming replay across all candidate rows. Numeric audit is bounded; Normal1 raw counts not tested here; biological donor assignment remains separate','source_identity_sha256':digest(matches)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--barcodes',type=Path,required=True)
    p.add_argument('--features',type=Path,required=True)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--preceding-barcodes',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():raise ContractError('refuse overwrite existing source audit')
    result=counts_fingerprint(a.matrix,a.barcodes,a.features,a.candidate,preceding_barcodes=a.preceding_barcodes)
    write_json(a.out,result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
if __name__=='__main__':
    main()
