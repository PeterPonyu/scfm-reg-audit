#!/usr/bin/env python
"""R2.2: where the "5.5 expected random motif hits per peak" comes from.

Scanner facts (src/v2/build_atac_graph_v2.py, src/v2/motif_utils.py): each linked
peak is resized to a 500 bp window around its midpoint; MOODS scans both strands
against a uniform background at p < 1e-5; scan_seq_hits returns the set of motifs
with >= 1 hit, so repeated or overlapping hits of one motif in a peak count once;
the peak x TF indicator then takes an OR over each TF's motifs, so a TF counts
once per peak however many of its motifs or positions hit. Heterodimer motifs
(A::B) count for both partners.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import log  # noqa: E402

sys.path.insert(0, str(common.SRC / "v2"))
import motif_utils as mu  # noqa: E402

ROOT = Path(common.fpa.ROOT)
L, P = 500, 1e-5
META = {"brain": "G_ATAC_v2_meta.json", "pbmc": "G_ATAC_v2_PBMC10k_meta.json",
        "fibroblast_mix": "G_ATAC_v2_GSE206767_meta.json"}


def main():
    mot = mu.parse_meme(str(ROOT / "data/motifs/JASPAR2024_CORE_vertebrates.meme"))
    man = json.loads(Path(common.fpa.MANI).read_text())
    panel = set(man["genes"])
    use, tf2mot = {}, {}
    for mid, d in mot.items():
        for t in mu.tf_symbols(d["name"]):
            if t in panel and man["tf_flags"].get(t):
                use[mid] = d
                tf2mot.setdefault(t, set()).add(mid)
    width = {mid: int(np.asarray(d["pwm"]).shape[1] if np.asarray(d["pwm"]).shape[0] == 4
                      else np.asarray(d["pwm"]).shape[0]) for mid, d in use.items()}
    w = np.array(list(width.values()))
    positions = 2 * (L - w + 1)                     # both strands
    first_order = len(use) * 2 * L * P
    e_positions = float(np.sum(positions * P))
    e_motifs = float(np.sum(1 - (1 - P) ** positions))
    e_tfs = float(sum(1 - np.prod([(1 - P) ** (2 * (L - width[m] + 1)) for m in ms])
                      for ms in tf2mot.values()))
    observed = {}
    for tissue, f in META.items():
        meta = json.loads((Path(common.fpa.OUT) / f).read_text())
        assert abs(meta["motif_p"] - P) < 1e-15, (tissue, meta["motif_p"])
        per_peak = meta["peak_motif_hits"] / meta["relevant_peaks"]
        observed[tissue] = {"linked_peaks": meta["relevant_peaks"],
                            "distinct_motif_hits": meta["peak_motif_hits"],
                            "per_peak": per_peak,
                            "ratio_to_expected": per_peak / e_motifs}
    k = [len(v) for v in tf2mot.values()]
    result = {
        "schema_version": 1,
        "analysis": "motif_hit_expectation",
        "jaspar_motifs_total": len(mot),
        "motifs_used": len(use),
        "panel_tfs_with_motif": len(tf2mot),
        "tfs_with_multiple_motifs": int(sum(x > 1 for x in k)),
        "max_motifs_per_tf": int(max(k)),
        "heterodimer_motifs_used": int(sum("::" in d["name"] for d in use.values())),
        "motif_width_bp": {"min": int(w.min()), "median": float(np.median(w)), "max": int(w.max())},
        "window_bp": L, "p_threshold": P, "background": "uniform (0.25 each base)",
        "expected_first_order": first_order,
        "expected_hit_positions": e_positions,
        "expected_distinct_motifs_per_peak": e_motifs,
        "expected_distinct_tfs_per_peak": e_tfs,
        "observed": observed,
        "counting_rule": ("distinct motifs per peak (any position, either strand); a TF counts "
                          "once per peak after an OR over its motifs; heterodimers count for both"),
    }
    common.write_json(common.OUT_DIR / "motif_expectation.json", result)
    log(f"motifs used {len(use)}/{len(mot)}; first-order {first_order:.2f}; "
        f"distinct motifs {e_motifs:.2f}; distinct TFs {e_tfs:.2f}")
    for t, o in observed.items():
        log(f"  {t}: {o['per_peak']:.2f}/peak = {o['ratio_to_expected']:.2f}x")


if __name__ == "__main__":
    main()
