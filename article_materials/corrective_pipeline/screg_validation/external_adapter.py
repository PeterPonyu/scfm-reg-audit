"""Score explicitly supplied assay labels against frozen per-donor model graphs."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .common import ContractError, cache_receipt, load_json, sha256
from .design import verify_panel
from .replication import source_signature
from .statistics import external_metrics, donor_interval, adjust, exact_paired_signflip, validate_external_labels


def score_labels(label_path, directory, readouts, *, expected_tissue, min_each, resamples, seed, expected_context=None):
    directory = Path(directory)
    panel = load_json(directory / "panel.json")
    prepared = load_json(directory / "prepared.json")
    verify_panel(panel)
    ctx = prepared["source_context"]
    if ctx.get("code") != source_signature() or (expected_context is not None and ctx != expected_context):
        raise ContractError("external scoring context drift")
    if not cache_receipt(directory / "prepared.receipt.json", ctx):
        raise ContractError("verified prepared receipt required for external scoring")
    positions = {g: i for i, g in enumerate(panel["genes"])}
    with Path(label_path).open() as f:
        labels = list(csv.DictReader(f))
    if not labels:
        raise ContractError("empty external assay label file")
    validate_external_labels(labels, expected_tissue=expected_tissue)
    report, summaries = [], []
    evaluation_donors = panel["split"]["evaluation"]
    for readout in readouts:
        scored, excluded = [], {"outside_panel": 0, "outside_tf_panel": 0, "non_evaluation_donor": 0, "self_edge": 0, "missing_model_score": 0}
        graphs = {}
        for donor in evaluation_donors:
            stem = prepared["donor_stems"][donor]
            task = directory / (stem + "-" + readout)
            if not cache_receipt(task / "result.receipt.json", {**ctx, "donor": donor, "readout": readout,
                                  "panel_sha256": panel["panel_sha256"]}):
                raise ContractError("verified model graph receipt required for external scoring")
            graphs[donor] = (np.load(directory / (stem + "-" + readout) / "model_graph.npz", allow_pickle=False)["graph"],
                             np.load(directory / (stem + "-graphs.npz"), allow_pickle=False)["coexpression"])
        for row in labels:
            if row.get("donor") not in evaluation_donors:
                excluded["non_evaluation_donor"] += 1
                continue
            if row.get("tf") not in positions or row.get("target") not in positions:
                excluded["outside_panel"] += 1
                continue
            if row["tf"] not in panel["tf_genes"]:
                excluded["outside_tf_panel"] += 1
                continue
            if row["tf"] == row["target"]:
                excluded["self_edge"] += 1
                continue
            i, j = positions[row["tf"]], positions[row["target"]]
            model, baseline = graphs[row["donor"]]
            if not np.isfinite(model[i, j]):
                excluded["missing_model_score"] += 1
                continue
            scored.append({**row, "model_score": float(model[i, j]), "baseline_score": float(baseline[i, j])})
        metrics = external_metrics(scored, expected_tissue=expected_tissue, min_each=min_each)
        report.append({"readout": readout, "metrics": metrics, "excluded": excluded})
        for task in ("verified_classification", "positive_vs_unlabeled_occupancy"):
            rows = [r for r in metrics if r["task"] == task]
            if not rows:
                continue
            for metric in ("AUC", "AP"):
                by_donor = {d: [r[metric]["paired_difference"] for r in rows if r["donor"] == d] for d in evaluation_donors}
                if any(not values for values in by_donor.values()):
                    summaries.append({"readout": readout, "task": task, "metric": metric, "donor_interval": None,
                                      "reason": "external label coverage absent in a frozen evaluation donor"})
                    continue
                differences = [float(np.mean(by_donor[d])) for d in evaluation_donors]
                summaries.append({"readout": readout, "task": task, "metric": metric,
                                  "donor_interval": donor_interval(differences, evaluation_donors,
                                         expected_donors=evaluation_donors, resamples=resamples, seed=seed),
                                  "paired_signflip_p": exact_paired_signflip(differences)})
    valid = [r for r in summaries if "paired_signflip_p" in r]
    if valid:
        p_family = [r["paired_signflip_p"] for r in valid] + [1.] * (len(readouts) * 4 - len(valid))
        for row, p in zip(valid, adjust(p_family, "Holm")[:len(valid)]):
            row["paired_signflip_Holm"] = float(p)
    return {"status": "complete", "readouts": report, "donor_summary": summaries,
            "label_sha256": sha256(label_path), "panel_sha256": panel["panel_sha256"],
            "registered_contrasts": len(readouts) * 4, "untestable_contrasts_retained_at_p1": len(readouts) * 4 - len(valid),
            "scope": "external assay transfer; true-negative and positive-vs-unlabeled tasks separated; assay labels do not validate necessity of five gates"}
