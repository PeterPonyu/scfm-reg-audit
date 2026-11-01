"""Comparable metrics, donor intervals and conditional component randomization."""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata, beta
from sklearn.metrics import roc_auc_score, average_precision_score

from .common import ContractError, identities, unique


def adjust(p, method="BY"):
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or not len(p) or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ContractError("invalid multiplicity family")
    order = np.argsort(p, kind="stable")
    m = len(p)
    if method == "BY":
        values = np.minimum.accumulate((p[order] * m * np.sum(1 / np.arange(1, m + 1)) / np.arange(1, m + 1))[::-1])[::-1]
    elif method == "Holm":
        values = np.maximum.accumulate(p[order] * np.arange(m, 0, -1))
    else:
        raise ContractError("unregistered correction")
    out = np.empty(m)
    out[order] = np.minimum(values, 1)
    return out


def partial_spearman(x, y, controls):
    if not np.isfinite(np.asarray(x, float)).all() or not np.isfinite(np.asarray(y, float)).all():
        raise ContractError("nonfinite partial-correlation scores")
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 1 or y.shape != x.shape or not len(x) or np.ptp(x) == 0 or np.ptp(y) == 0:
        raise ContractError("constant or invalid partial-correlation statistic")
    x, y = rankdata(x), rankdata(y)
    x, y = x - x.mean(), y - y.mean()
    controls = np.asarray(controls, float)
    if controls.ndim == 1:
        controls = controls[:, None]
    if x.ndim != 1 or y.shape != x.shape or controls.shape[0] != len(x) or not np.isfinite(controls).all():
        raise ContractError("partial correlation inputs misaligned")
    c = np.column_stack([np.ones(len(x)), controls])
    u, s, _ = np.linalg.svd(c, full_matrices=False)
    q = u[:, s > np.finfo(float).eps * max(c.shape) * s[0]]
    xr, yr = x - q @ (q.T @ x), y - q @ (q.T @ y)
    nx, ny = np.linalg.norm(xr), np.linalg.norm(yr)
    tolerance = np.finfo(float).eps * max(c.shape) * 8
    den = nx * ny
    if nx <= tolerance * np.linalg.norm(x) or ny <= tolerance * np.linalg.norm(y):
        raise ContractError("constant or fully conditioned statistic")
    return float(np.clip(xr @ yr / den, -1, 1))


def component_randomization(scores, graph, tf_indices, controls, *, strata, permutations, seed,
                            degree_mode="own", reference_degree=None):
    graph = np.asarray(graph, float)
    scores = np.asarray(scores, float)
    n = graph.shape[0]
    if graph.shape != (n, n) or scores.shape != graph.shape or not np.isfinite(graph).all() or not np.isfinite(scores).all():
        raise ContractError("invalid randomization graph")
    tf = np.asarray(tf_indices, int)
    if not len(tf) or len(tf) != len(set(tf)) or tf.min() < 0 or tf.max() >= n:
        raise ContractError("invalid TF rows")
    ii = np.repeat(tf, n)
    jj = np.tile(np.arange(n), len(tf))
    ii, jj = ii[ii != jj], jj[ii != jj]
    controls = np.asarray(controls, float)
    strata = np.asarray(strata)
    if controls.shape[0] != len(ii) or strata.shape != (n,) or permutations < 1:
        raise ContractError("randomization controls/strata misaligned")
    if degree_mode not in {"own", "reference", "none"}:
        raise ContractError("unknown degree estimand")
    if degree_mode == "reference" and (reference_degree is None or len(reference_degree) != 2):
        raise ContractError("independent degree source required")
    blocks = [np.flatnonzero(strata == label) for label in np.unique(strata)]
    if not any(len(b) > 1 for b in blocks):
        raise ContractError("all strata singleton; randomization uninformative")
    def statistic(g):
        c = controls
        if degree_mode != "none":
            od, ind = ((g > 0).sum(1), (g > 0).sum(0)) if degree_mode == "own" else reference_degree
            c = np.column_stack([c, np.asarray(od)[ii], np.asarray(ind)[jj]])
        return partial_spearman(scores[ii, jj], g[ii, jj], c)
    observed = statistic(graph)
    rng = np.random.default_rng(seed)
    nulls = {"M": [], "D": []}
    for _ in range(permutations):
        perm = np.arange(n)
        for b in blocks:
            perm[b] = rng.permutation(b)
        nulls["M"].append(statistic(graph[np.ix_(perm, perm)]))
        shuffled = graph.copy()
        for t in tf:
            for b in blocks:
                eligible = b[b != t]
                shuffled[t, eligible] = graph[t, rng.permutation(eligible)]
        nulls["D"].append(statistic(shuffled))
    p = {key: float((1 + (np.abs(values) >= abs(observed)).sum()) / (permutations + 1))
         for key, values in nulls.items()}
    return {"observed": observed, "pM": p["M"], "pD": p["D"], "pIUT": max(p.values()),
            "permutations": permutations, "degree_mode": degree_mode,
            "stratum_sizes": [len(b) for b in blocks], "nulls": nulls,
            "scope": "conditional reference distributions; real within-stratum exchangeability is not established by this test"}


def donor_interval(values, donor_ids, *, expected_donors, resamples=2000, seed=5613, min_donors=3):
    donors = unique(donor_ids, "evaluation donor")
    unique(expected_donors, "expected donor")
    values = np.asarray(values, float)
    if set(donors) != set(expected_donors) or values.shape != (len(donors),) or not np.isfinite(values).all():
        raise ContractError("donor interval requires the full frozen eligible donor set")
    point = float(values.mean())
    if len(donors) < min_donors:
        return {"point": point, "ci95": None, "n_donors": len(donors),
                "scope": "descriptive; too few independent donors for registered interval"}
    if resamples < 100:
        raise ContractError("at least 100 bootstrap draws required")
    rng = np.random.default_rng(seed)
    boot = values[rng.integers(0, len(values), size=(resamples, len(values)))].mean(1)
    return {"point": point, "ci95": list(map(float, np.quantile(boot, [.025, .975]))),
            "n_donors": len(donors), "resamples": resamples,
            "scope": "equal-donor percentile bootstrap; small donor counts limit interval precision"}


def exact_paired_signflip(values):
    values = np.asarray(values, float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all() or len(values) > 20:
        raise ContractError("exact sign-flip requires 1..20 finite donor differences")
    signs = 1 - 2 * ((np.arange(2 ** len(values))[:, None] >> np.arange(len(values))) & 1)
    null = (signs * values).mean(1)
    return float(np.mean(np.abs(null) >= abs(values.mean()) - 1e-14))


def validate_external_labels(rows, *, expected_tissue):
    identities_seen = set()
    for row in rows:
        key = tuple(identities([row.get(k) for k in ("donor", "tf", "target")], "external assay identity"))
        if key in identities_seen:
            raise ContractError("duplicate external donor/TF/target label")
        identities_seen.add(key)
        if row.get("tissue") != expected_tissue:
            raise ContractError("external reference is not matched to the registered tissue")
        if row.get("label_status") not in {"positive", "verified_negative", "unlabeled"} or not row.get("evidence_id"):
            raise ContractError("external label lacks assay status/evidence")


def external_metrics(rows, *, expected_tissue, min_each=5):
    """Unlabeled examples are never silently promoted to verified negatives."""
    validate_external_labels(rows, expected_tissue=expected_tissue)
    groups = {}
    for row in rows:
        key = (row["donor"], row["tf"], row["target"])
        for metric in ("model_score", "baseline_score"):
            if not np.isfinite(float(row[metric])):
                raise ContractError("nonfinite external score")
        groups.setdefault(key[:2], []).append(row)
    out = []
    for (donor, tf), items in sorted(groups.items()):
        status = {r["label_status"] for r in items}
        verified = [r for r in items if r["label_status"] != "unlabeled"]
        usable = verified if "verified_negative" in status else items
        y = np.array([int(r["label_status"] == "positive") for r in usable])
        if min(int(y.sum()), int((y == 0).sum())) < min_each:
            continue
        record = {"donor": donor, "tf": tf, "n_positive": int(y.sum()), "n_comparator": int((y == 0).sum()),
                  "task": "verified_classification" if "verified_negative" in status else "positive_vs_unlabeled_occupancy"}
        for name, scorer in (("AUC", roc_auc_score), ("AP", average_precision_score)):
            model = float(scorer(y, [r["model_score"] for r in usable]))
            baseline = float(scorer(y, [r["baseline_score"] for r in usable]))
            record[name] = {"model": model, "baseline": baseline, "paired_difference": model - baseline}
        out.append(record)
    if not out:
        raise ContractError("no external TF meets frozen positive/comparator coverage")
    return out


GATE_NAMES = ("dual_full", "concordance", "non_degree_sign", "readout_sign", "above_baseline")


def gate_validation(rows):
    """A model/readout-level causal success class is separate from edge occupancy."""
    if not rows:
        raise ContractError("no biological gate reference classes")
    keys = unique([r["readout_id"] for r in rows], "readout reference class")
    del keys
    out = []
    for row in rows:
        if row.get("reference_class") not in {"validated_success", "validated_failure"} or not row.get("class_evidence_id"):
            raise ContractError("gate construct validation requires independent validated readout classes")
        if set(row.get("gates", {})) != set(GATE_NAMES) or any(type(v) is not bool for v in row["gates"].values()):
            raise ContractError("five frozen boolean gates required")
        if row.get("class_assigned_without_gate_outcomes") is not True:
            raise ContractError("biological class was not assigned independently of gate outcomes")
        out.append(row)
    report = {}
    for gate in (*GATE_NAMES, "all_five"):
        positive = [r for r in out if r["reference_class"] == "validated_success"]
        negative = [r for r in out if r["reference_class"] == "validated_failure"]
        passes = lambda r: all(r["gates"].values()) if gate == "all_five" else r["gates"][gate]
        def binomial(k, n):
            return {"numerator": k, "denominator": n, "estimate": k / n if n else None,
                    "exact_ci95": [0. if k == 0 else float(beta.ppf(.025, k, n-k+1)),
                                   1. if k == n else float(beta.ppf(.975, k+1, n-k))] if n else None}
        report[gate] = {"sensitivity": binomial(sum(passes(r) for r in positive), len(positive)),
                        "specificity": binomial(sum(not passes(r) for r in negative), len(negative))}
    return {"metrics": report, "scope": "validation on independent readout reference classes; correlation does not prove biological necessity"}
