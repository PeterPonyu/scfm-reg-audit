#!/usr/bin/env python
"""Split the combined 10x cellranger-arc PBMC Multiome h5 (RNA+ATAC, SAME cells) into
RNA h5ad (var_names=gene symbol) and ATAC h5ad (var_names='chr:start-end'), same barcodes —
matching the formats build_atac_graph_v2.py / fm_readout.py already expect."""
import scipy.io, scipy.sparse as sp, numpy as np, anndata as ad, h5py, pandas as pd, time
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
SRC = f"{ROOT}/data/multiome/pbmc10k_multiome_filtered.h5"
f = h5py.File(SRC, "r"); m = f["matrix"]
shape = m["shape"][:]  # (n_features, n_cells)
X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(shape)).T.tocsr()  # cells x features
bc = np.array([b.decode() for b in m["barcodes"][:]])
ft = np.array([b.decode() for b in m["features"]["feature_type"][:]])
name = np.array([b.decode() for b in m["features"]["name"][:]])
fid = np.array([b.decode() for b in m["features"]["id"][:]])
log(f"loaded combined matrix {X.shape} (cells x features), RNA={int((ft=='Gene Expression').sum())} ATAC={int((ft=='Peaks').sum())}")

rna_i = np.where(ft == "Gene Expression")[0]
atac_i = np.where(ft == "Peaks")[0]

Arna = ad.AnnData(X=X[:, rna_i].tocsr(), obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=name[rna_i]))
Arna.var_names_make_unique()
Arna.write_h5ad(f"{ROOT}/data/multiome/pbmc10k_rna.h5ad")
log(f"SAVED pbmc10k_rna.h5ad {Arna.shape}")

Aatac = ad.AnnData(X=X[:, atac_i].tocsr(), obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=fid[atac_i]))
Aatac.write_h5ad(f"{ROOT}/data/multiome/pbmc10k_atac.h5ad")
log(f"SAVED pbmc10k_atac.h5ad {Aatac.shape}")
