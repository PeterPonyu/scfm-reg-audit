#!/usr/bin/env python
"""Build pair-level (TF, gene) probe features from pooled graph caches.

Each sample is an ordered pair (tf, gene). A family's feature block is derived
from that family's gene x gene graph only, so the probe generalises to TFs it
never saw: nothing in the design is indexed by TF identity.

Feature block per family (5 columns):
  0  w        edge weight  S[tf, gene]
  1  w_rev    reverse edge S[gene, tf]  (0 if the graph is symmetric)
  2  tf_deg   row mean of the TF, broadcast over its genes
  3  g_deg    column mean of the target gene
  4  rank_g   within-TF rank of w among the 1200 genes, scaled to [0, 1]

Targets are the pooled PBMC ATAC-proxy rows.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPLIT = ROOT / "results/v2/tf_disjoint_split_v2.json"
G_ATAC = ROOT / "results/v2/G_ATAC_v2_PBMC10k.npz"
MANIFEST = ROOT / "data/manifest/shared_genes.v2.json"
CONFOUNDS = ROOT / "results/v2/pbmc_confounds_v2.npz"

PBMC_POOLED = ROOT / "results/v2/pbmc_fmgraphs_pooled.npz"
PBMC_SCGPT = ROOT / "results/v2/pbmc_scgpt_pooled_v2.npz"
PBMC_UCE = ROOT / "results/v2/pbmc_uce_pooled_v2.npz"
BRAIN_FLOOR = ROOT / "results/v2/brain_floor_graph_v2.npz"

OUT_DIR = ROOT / "results/v2/tf_probe_pair"
N_FEAT = 5


def load_graphs(genes):
    """Return {family: gene x gene matrix}, all on the manifest gene order."""
    graphs = {}
    with np.load(PBMC_POOLED, allow_pickle=False) as d:
        graphs["co_expression"] = np.asarray(d["co"], dtype=np.float64)
        graphs["geneformer_embed"] = np.asarray(d["gf"], dtype=np.float64)
        graphs["geneformer_attn"] = np.asarray(d["at"], dtype=np.float64)
    with np.load(PBMC_SCGPT, allow_pickle=False) as d:
        if [str(g) for g in d["genes"]] != genes:
            raise ValueError("scGPT cache gene order does not match manifest")
        graphs["scGPT_encoder"] = np.asarray(d["sg"], dtype=np.float64)
    with np.load(PBMC_UCE, allow_pickle=False) as d:
        if [str(g) for g in d["genes"]] != genes:
            raise ValueError("UCE cache gene order does not match manifest")
        graphs["UCE_encoder"] = np.asarray(d["uce"], dtype=np.float64)
    with np.load(BRAIN_FLOOR, allow_pickle=False) as d:
        graphs["random_floor"] = np.asarray(d["G"], dtype=np.float64)

    n = len(genes)
    for label, mat in graphs.items():
        if mat.shape != (n, n):
            raise ValueError(f"{label}: shape {mat.shape} != ({n}, {n})")
    return graphs


def pair_block(mat, tf_rows):
    """Featurise (tf, gene) pairs for one graph. Returns (n_tf * n_gene, N_FEAT)."""
    sub = mat[tf_rows, :]                      # (n_tf, n_gene)
    rev = mat[:, tf_rows].T                    # (n_tf, n_gene)
    tf_deg = np.repeat(sub.mean(axis=1, keepdims=True), sub.shape[1], axis=1)
    g_deg = np.repeat(mat.mean(axis=0, keepdims=True), sub.shape[0], axis=0)
    order = np.argsort(np.argsort(sub, axis=1), axis=1)
    rank_g = order / max(sub.shape[1] - 1, 1)
    return np.stack([sub, rev, tf_deg, g_deg, rank_g], axis=-1).reshape(-1, N_FEAT)


def main():
    """Write per-family pair-level feature blocks plus targets and confounds."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    genes = json.loads(MANIFEST.read_text())["genes"]
    n_genes = len(genes)
    split = json.loads(SPLIT.read_text())
    train_idx = np.asarray(split["train_idx"], dtype=int)
    test_idx = np.asarray(split["test_idx"], dtype=int)
    if set(train_idx.tolist()) & set(test_idx.tolist()):
        raise ValueError("train/test TF indices overlap")

    with np.load(G_ATAC, allow_pickle=False) as d:
        if [str(g) for g in d["genes"]] != genes:
            raise ValueError("G_ATAC gene order does not match manifest")
        types = [str(t) for t in d["types"]]
        G = np.mean([np.asarray(d[f"G_{t}"], dtype=np.float64) for t in types], axis=0)
        tf_rows_all = d["tf_rows"].astype(int)

    tf_rows_train = tf_rows_all[train_idx]
    tf_rows_test = tf_rows_all[test_idx]
    y_train = G[tf_rows_train, :].reshape(-1)
    y_test = G[tf_rows_test, :].reshape(-1)

    with np.load(CONFOUNDS, allow_pickle=False) as d:
        if [str(g) for g in d["genes"]] != genes:
            raise ValueError("confound gene order does not match manifest")
        conf_gene = np.column_stack(
            [np.asarray(d[k], dtype=np.float64) for k in ("peakcount", "genelen", "gc", "detv")]
        )

    graphs = load_graphs(genes)
    print(f"{n_genes} genes | train TFs {len(train_idx)} | test TFs {len(test_idx)}")

    for label, mat in graphs.items():
        xtr = pair_block(mat, tf_rows_train)
        xte = pair_block(mat, tf_rows_test)
        if xtr.shape != (len(train_idx) * n_genes, N_FEAT):
            raise ValueError(f"{label}: train block shape {xtr.shape}")
        np.savez_compressed(
            OUT_DIR / f"{label}_pairs.npz",
            X_train=xtr.astype(np.float32),
            X_test=xte.astype(np.float32),
            feature_names=["w", "w_rev", "tf_deg", "g_deg", "rank_g"],
            family=label,
        )
        print(f"  {label:20s} train {xtr.shape} test {xte.shape}")

    np.savez_compressed(
        OUT_DIR / "pair_targets.npz",
        y_train=y_train.astype(np.float32),
        y_test=y_test.astype(np.float32),
        train_idx=train_idx,
        test_idx=test_idx,
        train_tfs=np.asarray(split["train_tfs"], dtype="<U32"),
        test_tfs=np.asarray(split["test_tfs"], dtype="<U32"),
        genes=np.asarray(genes, dtype="<U32"),
        conf_gene=conf_gene.astype(np.float32),
        n_genes=n_genes,
    )
    prov = {
        "schema_version": 1,
        "design": "pair-level (TF, gene) probe; TF-disjoint train/test",
        "families": sorted(graphs),
        "feature_names": ["w", "w_rev", "tf_deg", "g_deg", "rank_g"],
        "n_genes": n_genes,
        "n_train_tfs": int(len(train_idx)),
        "n_test_tfs": int(len(test_idx)),
        "n_train_pairs": int(len(y_train)),
        "n_test_pairs": int(len(y_test)),
        "split_seed": split["split_seed"],
        "atac_celltypes": types,
    }
    (OUT_DIR / "provenance.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    print(f"targets: train {y_train.shape} test {y_test.shape}")


if __name__ == "__main__":
    main()
