#!/usr/bin/env python
"""Split a 10x cellranger-arc H5 into RNA/ATAC H5AD with the same barcodes."""
import argparse
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))


def split_multiome(source, out_dir):
    source, out_dir = Path(source), Path(out_dir)
    if not source.is_file():
        raise FileNotFoundError(f"10x multiome input missing: {source}")
    with h5py.File(source, "r") as handle:
        m = handle["matrix"]
        shape = tuple(m["shape"][:])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=shape).T.tocsr()
        bc = np.array([b.decode() for b in m["barcodes"][:]])
        ft = np.array([b.decode() for b in m["features"]["feature_type"][:]])
        name = np.array([b.decode() for b in m["features"]["name"][:]])
        fid = np.array([b.decode() for b in m["features"]["id"][:]])
    rna_i, atac_i = np.where(ft == "Gene Expression")[0], np.where(ft == "Peaks")[0]
    if not len(rna_i) or not len(atac_i):
        raise ValueError("Expected both Gene Expression and Peaks features")
    rna = ad.AnnData(X=X[:, rna_i], obs=pd.DataFrame(index=bc),
                    var=pd.DataFrame(index=name[rna_i]))
    rna.var_names_make_unique()
    atac = ad.AnnData(X=X[:, atac_i], obs=pd.DataFrame(index=bc),
                     var=pd.DataFrame(index=fid[atac_i]))
    out_dir.mkdir(parents=True, exist_ok=True)
    rna_path, atac_path = out_dir / "pbmc10k_rna.h5ad", out_dir / "pbmc10k_atac.h5ad"
    rna.write_h5ad(rna_path)
    atac.write_h5ad(atac_path)
    return rna_path, atac_path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=ROOT / "data/multiome/pbmc10k_multiome_filtered.h5")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "data/multiome")
    args = ap.parse_args()
    for path in split_multiome(args.input, args.out_dir):
        print(f"SAVED {path}")


if __name__ == "__main__":
    main()
