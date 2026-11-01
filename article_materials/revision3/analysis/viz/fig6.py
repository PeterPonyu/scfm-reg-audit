"""Figure 6: donor-separated bone-marrow extension and external occupancy references."""
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from sklearn.metrics import average_precision_score, roc_auc_score

from . import style as S
from .common import (FIGDATA, HIST, READOUTS, ROOT, clean, fmt, group_labels, legend_above,
                     load, payload, pt, short)
from .layout import Canvas, split


class Data:
    def __init__(self):
        self.bmmc = payload("figure5_independent_validation")
        self.rows = self.bmmc["rows"]
        self.final = self.bmmc["final"]["rows"]
        self.prepared = self.bmmc["prepared"]
        self.donors = self.prepared["split"]["evaluation"]
        pool = load(HIST / "validation/qa/team_20260930/red_donor_coverage_potential_fullpool.json")["donors"]
        self.pool = [next(r for r in pool if r["donor"] == dn) for dn in self.donors]
        self.ctcf = load(ROOT / "inputs/ctcf-compact.json")
        self.external = load(FIGDATA / "external_occupancy.json")
        self.ctcf_metrics()

    def ctcf_metrics(self):
        y = np.array([v == "positive" for v in self.ctcf["labels"]], int)
        scores = {"co_expression": self.ctcf["baseline_scores"],
                  **{r["readout"]: r["scores"] for r in self.ctcf["readouts"] if r["external_eligible"]}}
        metrics = []
        for name, values in scores.items():
            v = np.array(values, float)
            if np.isfinite(v).all():
                metrics.append({"readout": name, "auc": float(roc_auc_score(y, v)), "ap": float(average_precision_score(y, v))})
        self.y = y
        self.base = next(m for m in metrics if m["readout"] == "co_expression")
        self.models = [m for m in metrics if m is not self.base]
        assert len(y) == 1199 and int(y.sum()) == 1016 and len(self.models) == 5


def panel_a(d):
    def draw(ax):
        split_ = d.prepared["split"]
        roles = [("Discovery", S.BLUE, split_["discovery"]), ("Degree reference", S.ORANGE, split_["degree_reference"]),
                 ("Evaluation", S.TEAL, split_["evaluation"])]
        assert [len(r[2]) for r in roles] == [3, 2, 5]
        x, rec = 0.0, {}
        for name, color, donors in roles:
            start = x
            for dn in donors:
                ax.add_patch(Rectangle((x, 0), 0.94, 0.9, facecolor=color, edgecolor="none"))
                ax.text(x + 0.47, 0.45, str(dn), ha="center", va="center", color="white", fontsize=S.FS_NOTE)
                x += 1
            ax.text((start + x - 0.06) / 2, 1.18, name, ha="center", va="bottom", fontsize=S.FS_NOTE, color=color, fontweight="normal")
            rec[name] = list(donors)
            x += 0.35
        ax.set(xlim=(0, x - 0.35), ylim=(0, 1.6))
        ax.axis("off")
        return {"roles": rec}
    return draw


def panel_b(d):
    def draw(ax):
        mat = np.array([[next(r["pair_coverage"] for r in d.rows if r["donor"] == dn and r["readout"] == m) for m in READOUTS]
                        for dn in d.donors])
        complete = mat == 1
        ax.imshow(complete, cmap=ListedColormap([S.LIGHT, S.TEAL]), vmin=0, vmax=1, aspect="auto")
        for i in range(5):
            for j in range(6):
                ax.text(j, i, f"{100 * mat[i, j]:.1f}", ha="center", va="center", fontsize=S.FS_NOTE,
                        color="white" if complete[i, j] else "black")
        for k in range(1, 5):
            ax.axhline(k - 0.5, color="white", lw=0.7)
        for k in range(1, 6):
            ax.axvline(k - 0.5, color="white", lw=0.7)
        ax.set_yticks(range(5)); ax.set_yticklabels(d.donors)
        ax.set_xticks(range(6)); ax.set_xticklabels([short(m) for m in READOUTS], rotation=90)
        ax.set_ylabel("Donor")
        ax.tick_params(length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        return {"donors": d.donors, "readouts": [short(m) for m in READOUTS], "pair_coverage_pct": (100 * mat).tolist()}
    return draw


def panel_c(d):
    def draw(ax):
        counts = [sum(r["eligible"] for r in d.rows if r["donor"] == dn) for dn in d.donors]
        ax.bar(range(5), counts, color=S.TEAL, width=0.6)
        ax.bar(range(5), 6 - np.array(counts), bottom=counts, color=S.LIGHT, edgecolor=S.GRAY, lw=0.5, width=0.6)
        for i, v in enumerate(counts):
            ax.text(i, 6.3, f"{v}/6", ha="center", fontsize=S.FS_NOTE)
        ax.set(xlim=(-0.6, 4.6), ylim=(0, 7.7), yticks=[0, 3, 6], ylabel="Admitted tasks")
        ax.set_xticks(range(5)); ax.set_xticklabels(d.donors, rotation=60, ha="right", rotation_mode="anchor")
        clean(ax, grid="y")
        return {"donors": d.donors, "admitted": counts, "tasks_per_donor": 6}
    return draw


def panel_d(d):
    def draw(ax):
        names = [m for m in READOUTS if any(r["readout"] == m for r in d.final)]
        assert len(names) == 5
        rec = {}
        for deg, color, marker, off in (("reference", S.BLUE, "o", -0.25), ("own", S.TEAL, "s", 0.0), ("none", S.ORANGE, "^", 0.25)):
            vals = [next(r["observed"] for r in d.final if r["readout"] == m and r["degree_mode"] == deg and r["null_mode"] == "restricted")
                    for m in names]
            ax.scatter(vals, np.arange(5) + off, color=color, marker=marker, s=19, zorder=3, lw=0)
            rec[deg] = vals
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(ylim=(4.6, -0.6), xlim=(-0.03, 0.023), xticks=[-0.02, 0, 0.02], xlabel="Partial Spearman ρ")
        ax.set_xticklabels(["−0.02", "0", "0.02"])
        ax.set_yticks(range(5)); ax.set_yticklabels([short(m) for m in names])
        clean(ax)
        legend_above(ax, [pt("o", S.BLUE, "Ref."), pt("s", S.TEAL, "Own"), pt("^", S.ORANGE, "None")])
        return {"readouts": [short(m) for m in names], "observed": rec, "donor": "16710"}
    return draw


def panel_e(d):
    def draw(ax):
        missing = [r["attention_missing_pairs"] for r in d.pool]
        ys = np.arange(len(missing))
        ax.barh(ys, missing, height=0.6, color=S.ORANGE)
        for y, v in zip(ys, missing):
            ax.text(v + 600, y, f"{v:,}", va="center", fontsize=S.FS_NOTE)
        cells = [f"{r['cells']:,}" for r in d.pool]
        absent = [str(len(r["full_pool_absent_TFs"])) for r in d.pool]
        for y, c_, a_ in zip(ys, cells, absent):
            ax.text(0.80, y, c_, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=S.FS_NOTE)
            ax.text(0.965, y, a_, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=S.FS_NOTE)
        ax.text(0.80, -0.95, "Cells", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=S.FS_NOTE, fontweight="normal")
        ax.text(0.965, -0.95, "Absent TFs", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=S.FS_NOTE, fontweight="normal")
        ax.set(xlim=(0, 55000), ylim=(len(ys) - 0.45, -1.3), xticks=[0, 10000, 20000, 30000],
               xlabel="Attention pairs unavailable from any subset of the pool")
        ax.set_xticklabels(["0", "10k", "20k", "30k"])
        ax.set_yticks(ys); ax.set_yticklabels(d.donors)
        clean(ax)
        return {"donors": d.donors, "attention_pairs_missing": missing, "cells": [r["cells"] for r in d.pool],
                "absent_tfs": [len(r["full_pool_absent_TFs"]) for r in d.pool]}
    return draw


def panel_f(d):
    def draw(ax):
        pos, n = int(d.y.sum()), len(d.y)
        ax.barh(0, pos, height=0.7, color=S.TEAL)
        ax.barh(0, n - pos, left=pos, height=0.7, color=S.LIGHT, edgecolor=S.GRAY, lw=0.5)
        ax.text(pos / 2, 0, f"{pos:,} supported by both ENCODE experiments", ha="center", va="center", color="white", fontsize=S.FS_NOTE)
        ax.text(pos + (n - pos) / 2, 0, f"{n - pos} unlabelled", ha="center", va="center", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, n), ylim=(-0.5, 0.5))
        ax.axis("off")
        return {"targets": n, "supported": pos, "unlabelled": n - pos}
    return draw


def metric_panel(d, key, label_x, xlim, xticks, show_labels):
    def draw(ax):
        ordered = d.models + [d.base]
        for i, r in enumerate(ordered):
            color = S.GRAY if r is d.base else S.TEAL
            ax.scatter(r[key], i, c=color, s=20, zorder=3, lw=0)
            delta = "" if r is d.base else f" ({fmt(r[key] - d.base[key], 3, sign=True)})"
            ax.text(r[key] + (xlim[1] - xlim[0]) * 0.035, i, f"{r[key]:.3f}{delta}", va="center", fontsize=S.FS_NOTE)
        ax.axvline(d.base[key], color=S.GRAY, ls="--", lw=0.7)
        ax.set(ylim=(5.6, -0.6), xlim=xlim, xticks=xticks, xlabel="ROC-AUC" if key == "auc" else "Average precision")
        ax.set_yticks(range(6))
        ax.set_yticklabels([short(r["readout"]) if r is not d.base else "Co-expr." for r in ordered] if show_labels else [])
        if not show_labels:
            ax.tick_params(axis="y", length=0)
        clean(ax)
        return {"readouts": [short(r["readout"]) if r is not d.base else "Co-expression" for r in ordered],
                key: [r[key] for r in ordered]}
    return draw


def historical(d):
    return [(t, r) for t in ("brain", "pbmc") for r in d.external["tissues"][t]["rows"] if "auc_delta_vs_co" in r]


def interval_panel(d, key, xlim, xticks, label_x, show_labels):
    def draw(ax):
        rows = historical(d)
        assert len(rows) == 13
        for i, (t, r) in enumerate(rows):
            color = S.BLUE if t == "brain" else S.TEAL
            v = r[key + "_delta_vs_co"]
            lo, hi = r[key + "_delta_tf_resampling95"]
            ax.errorbar(v, i, xerr=[[v - lo], [hi - v]], fmt="o", color=color, ms=3.2, capsize=1.8, lw=0.8)
            ax.text(label_x, i, f"{r['mean_' + key]:.3f}", va="center", fontsize=S.FS_NOTE)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.axhline(7.5, color=S.RULE, lw=0.6)
        labels = [short(r["readout"]) for _, r in rows]
        ax.set(ylim=(12.6, -0.6), xlim=xlim, xticks=xticks,
               xlabel="ΔAUC; right: mean" if key == "auc" else "ΔAP; right: mean")
        ax.set_xticklabels([fmt(v, 2 if key == "auc" else 2) if v else "0" for v in xticks])
        ax.set_yticks(range(13)); ax.set_yticklabels(labels if show_labels else [])
        if not show_labels:
            ax.tick_params(axis="y", length=0)
        clean(ax)
        if show_labels:
            group_labels(ax, [("Mixed", 0, 7), ("PBMC", 8, 12)], labels)
        return {"rows": [{"context": "Mixed" if t == "brain" else "PBMC", "readout": short(r["readout"]),
                          "delta": r[key + "_delta_vs_co"], "interval95": r[key + "_delta_tf_resampling95"],
                          "mean": r["mean_" + key]} for t, r in rows], "n_tfs": d.external["tissues"]["brain"]["n_tfs"]}
    return draw


def build(out_height=156.0):
    d = Data()
    c = Canvas("Figure6", out_height)
    ra, r1, r2, rf, r3 = split(c.rect, "v", [11.0, 37.0, 25.0, 9.0, 60.0])
    b_cell, c_cell, d_cell = split(r1, "h", [1.25, 0.8, 1.05])
    gh, ij = split(r3, "h", [1.12, 1.68])
    g_cell, h_cell = split(gh, "v", [1.0, 1.0])
    i_cell, j_cell = split(ij, "h", [1.2, 1.0])
    c.panel(ra, "A", "Frozen donor roles", panel_a(d), header=True)
    c.panel(b_cell, "B", "Pair coverage (%)", panel_b(d))
    c.panel(c_cell, "C", "Admitted tasks", panel_c(d))
    c.panel(d_cell, "D", "Donor 16710 effects", panel_d(d))
    c.panel(r2, "E", "Full assigned-pool raw-support ceiling", panel_e(d))
    c.panel(rf, "F", "CTCF targets", panel_f(d))
    c.panel(g_cell, "G", "CTCF ROC-AUC", metric_panel(d, "auc", 0, (0.38, 0.72), [0.4, 0.5, 0.6], True), align_left="gh")
    c.panel(h_cell, "H", "CTCF average precision", metric_panel(d, "ap", 0, (0.79, 0.955), [0.8, 0.85, 0.9], True), align_left="gh")
    c.panel(i_cell, "I", "204-TF ΔAUC", interval_panel(d, "auc", (-0.115, 0.075), [-0.1, 0], 0.012, True))
    c.panel(j_cell, "J", "204-TF ΔAP", interval_panel(d, "ap", (-0.05, 0.034), [-0.04, 0], 0.006, False))
    return c, {}
