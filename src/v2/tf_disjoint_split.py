#!/usr/bin/env python
"""Generate stratified train/test split for TF-disjoint probe experiment."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
G_ATAC = ROOT / "results/v2/G_ATAC_v2_PBMC10k.npz"
MANIFEST = ROOT / "data/manifest/shared_genes.v2.json"
OUT = ROOT / "results/v2/tf_disjoint_split_v2.json"

SPLIT_SEED = 20260730
TRAIN_FRAC = 0.7


def main():
    """Stratify TFs by outdegree quantiles, split 70/30 train/test."""
    with np.load(G_ATAC, allow_pickle=False) as data:
        genes = [str(g) for g in data["genes"]]
        tf_row_indices = data["tf_rows"].astype(int)  # Already indices
        types = [str(t) for t in data["types"]]
        # Pool across cell types: stack all G matrices and average
        G_list = [data[f"G_{ct}"] for ct in types]
        G_pooled = np.mean(G_list, axis=0)

    manifest = json.loads(MANIFEST.read_text())
    panel_genes = manifest["genes"]
    if genes != panel_genes:
        raise ValueError("G_ATAC gene order does not match manifest")

    n_tfs = len(tf_row_indices)
    n_train = int(n_tfs * TRAIN_FRAC)

    # Extract TF gene names from indices
    tf_names = [genes[i] for i in tf_row_indices]

    # Compute TF outdegree: for each TF row, count how many genes it regulates
    G_tf = G_pooled[tf_row_indices, :]
    outdeg = np.sum(G_tf > 0, axis=1)
    if outdeg.shape[0] != n_tfs:
        raise ValueError(f"TF outdegree shape mismatch: {outdeg.shape[0]} vs {n_tfs} TFs")

    # Stratify by outdegree quartiles
    quartiles = np.percentile(outdeg, [25, 50, 75])
    strata = np.digitize(outdeg, quartiles)

    rng = np.random.RandomState(SPLIT_SEED)
    train_idx = []
    test_idx = []

    for stratum in np.unique(strata):
        stratum_mask = strata == stratum
        stratum_tfs = np.where(stratum_mask)[0]
        rng.shuffle(stratum_tfs)

        n_stratum_train = int(len(stratum_tfs) * TRAIN_FRAC)
        train_idx.extend(stratum_tfs[:n_stratum_train].tolist())
        test_idx.extend(stratum_tfs[n_stratum_train:].tolist())

    train_idx = sorted(train_idx)
    test_idx = sorted(test_idx)

    # Verify disjointness
    if set(train_idx) & set(test_idx):
        raise ValueError("Train/test overlap detected")
    if len(train_idx) + len(test_idx) != n_tfs:
        raise ValueError(f"Split size mismatch: {len(train_idx)} + {len(test_idx)} != {n_tfs}")

    train_tfs = [tf_names[i] for i in train_idx]
    test_tfs = [tf_names[i] for i in test_idx]
    train_outdeg = outdeg[train_idx]
    test_outdeg = outdeg[test_idx]

    split = {
        "schema_version": 1,
        "split_seed": SPLIT_SEED,
        "train_frac": TRAIN_FRAC,
        "n_tfs": n_tfs,
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "train_tfs": train_tfs,
        "test_tfs": test_tfs,
        "train_idx": train_idx,
        "test_idx": test_idx,
        "train_outdeg_mean": float(train_outdeg.mean()),
        "train_outdeg_std": float(train_outdeg.std()),
        "test_outdeg_mean": float(test_outdeg.mean()),
        "test_outdeg_std": float(test_outdeg.std()),
        "quartiles": quartiles.tolist(),
        "provenance": {
            "G_ATAC": str(G_ATAC.relative_to(ROOT)),
            "manifest": str(MANIFEST.relative_to(ROOT)),
            "G_pooled_shape": list(G_pooled.shape),
            "n_cell_types": len(types),
        }
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(split, indent=2, sort_keys=True) + "\n")

    print(f"TF split: {len(train_idx)} train, {len(test_idx)} test")
    print(f"Train outdegree: {train_outdeg.mean():.1f} ± {train_outdeg.std():.1f}")
    print(f"Test outdegree: {test_outdeg.mean():.1f} ± {test_outdeg.std():.1f}")
    print(f"Saved to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
