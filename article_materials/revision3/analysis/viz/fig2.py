"""Figure 2: genuine-brain effects, permutation resolution and training-exposure limits."""
import numpy as np
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from . import style as S
from .common import (ADMITTED, DEGREE, FIGDATA, ROOT, clean, corrected_rows, fmt, hbars, heat, legend_above,
                     load, payload, pt, short)
from .layout import Canvas, split

NULLS = ("restricted", "unrestricted")


class Data:
    def __init__(self):
        self.rows = corrected_rows()
        self.resolution = load(FIGDATA / "resolution_n9999.json")
        self.training = payload("figure7_training")
        self.levels = load(ROOT / "inputs/exposure-evidence-levels.json")
        assert len(self.rows) == 16

    def pick(self, readout, degree, null):
        return next(r for r in self.rows if (r["readout"], r["degree_mode"], r["null_mode"]) == (readout, degree, null))

    def ordered(self):
        return [self.pick(m, deg, null) for deg in ("none", "own") for null in NULLS for m in ADMITTED]


def _dedupe_observed(d):
    """Restricted and unrestricted tests share one observed effect; assert before plotting it once."""
    for m in ADMITTED:
        for deg in ("none", "own"):
            a, b = d.pick(m, deg, "restricted"), d.pick(m, deg, "unrestricted")
            assert a["observed"] == b["observed"] and a["comparable_difference"] == b["comparable_difference"]


def panel_a(d):
    def draw(ax):
        _dedupe_observed(d)
        rec = {}
        for deg, color, marker, name in DEGREE:
            off = -0.25 if deg == "none" else 0.25
            vals = [d.pick(m, deg, "restricted")["observed"] for m in ADMITTED]
            ax.scatter(vals, np.arange(4) + off, c=color, marker=marker, s=24, zorder=3)
            rec[deg] = vals
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(yticks=range(4), ylim=(3.6, -0.6), xlim=(-0.0092, 0.0046), xticks=[-0.008, -0.004, 0, 0.004],
               xlabel="Partial Spearman ρ")
        ax.set_yticklabels([short(m) for m in ADMITTED])
        ax.set_xticklabels(["−0.008", "−0.004", "0", "0.004"])
        clean(ax)
        legend_above(ax, [pt("^", S.ORANGE, "No degree"), pt("s", S.TEAL, "Own degree")])
        return {"readouts": [short(m) for m in ADMITTED], "observed": rec}
    return draw


def panel_b(d):
    def draw(ax):
        rec = []
        for i, m in enumerate(ADMITTED):
            for deg, color, marker, name in DEGREE:
                r = d.pick(m, deg, "restricted")
                y = i + (-0.25 if deg == "none" else 0.25)
                ax.plot([r["comparable_fm"], r["comparable_baseline"]], [y, y], color=color, lw=0.9)
                ax.scatter(r["comparable_fm"], y, color=color, marker=marker, s=22, zorder=3)
                ax.scatter(r["comparable_baseline"], y, facecolor="white", edgecolor=color, s=24, zorder=3, lw=0.9)
                ax.text(0.0555, y, fmt(r["comparable_difference"], 4, sign=True), va="center", fontsize=S.FS_NOTE)
                rec.append({"readout": short(m), "degree": deg, "model": r["comparable_fm"],
                            "coexpression": r["comparable_baseline"], "difference": r["comparable_difference"]})
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(yticks=range(4), ylim=(3.6, -0.6), xlim=(-0.016, 0.082), xticks=[0, 0.02, 0.04],
               xlabel="Partial ρ (label: model − co-expr.)")
        ax.tick_params(axis="y", labelleft=False)
        clean(ax)
        legend_above(ax, [pt("o", S.GRAY, "Model"), pt("o", S.GRAY, "Co-expression", hollow=True)])
        return {"pairs": rec}
    return draw


def _slot_labels(rows):
    return [f"{short(r['readout'])} {'N' if r['degree_mode'] == 'none' else 'O'}-{r['null_mode'][0].upper()}"
            for r in rows]


def panel_c(d):
    def draw(ax):
        order = d.ordered()
        matrix = np.array([[r["pM"], r["pD"], r["pIUT"]] for r in order])
        heat(ax, matrix, cmap="Blues", norm=LogNorm(vmin=0.01, vmax=1), fmt="{:.4f}", dark=lambda v: v > 0.25,
             row_labels=_slot_labels(order), col_labels=["pM", "pD", "pIUT"], separators=(4, 8, 12))
        return {"slots": _slot_labels(order), "columns": ["pM", "pD", "pIUT"], "p_9999": matrix.tolist()}
    return draw


def panel_d(d):
    def draw(ax):
        order = d.ordered()
        h = 0.05 / 72
        by = h / sum(1 / np.arange(1, 73))
        for i, r in enumerate(order):
            ax.plot([r["registered_999"]["pIUT"], r["pIUT"]], [i, i], color=S.GRAY, lw=0.8)
            ax.scatter(r["registered_999"]["pIUT"], i - 0.12, facecolor="white", edgecolor=S.GRAY, s=18, zorder=3, lw=0.8)
            ax.scatter(r["pIUT"], i + 0.12, color=S.BLUE, s=18, zorder=3)
        ax.axvline(h, color=S.RED, ls="--", lw=0.8)
        ax.axvline(by, color=S.ORANGE, ls=":", lw=1.0)
        for y in (3.5, 7.5, 11.5):
            ax.axhline(y, color=S.RULE, lw=0.6)
        ax.set(xscale="log", xlim=(0.0001, 1.3), ylim=(15.5, -0.5), yticks=range(16), xticks=[1e-4, 1e-3, 1e-2, 1e-1, 1],
               xlabel="Raw IUT probability")
        ax.set_xticklabels(["0.0001", "0.001", "0.01", "0.1", "1"])
        ax.tick_params(axis="y", labelleft=False)
        ax.minorticks_off(); clean(ax)
        legend_above(ax, [pt("o", S.GRAY, "999", hollow=True), pt("o", S.BLUE, "9,999"),
                          Line2D([], [], color=S.RED, ls="--", label="Holm"), Line2D([], [], color=S.ORANGE, ls=":", label="BY")])
        return {"registered_999": [r["registered_999"]["pIUT"] for r in order], "post_hoc_9999": [r["pIUT"] for r in order],
                "holm_first_step": h, "by_first_step": by, "family": 72, "measured": 16}
    return draw


def panel_e(d):
    def draw(ax):
        rr = d.resolution["fm_rows"]
        x = [r["frozen_p"] for r in rr]
        y = [r["p"] for r in rr]
        ax.plot([8e-5, 1.1], [8e-5, 1.1], color=S.GRAY, lw=0.7)
        ax.scatter(x, y, s=14, color=S.BLUE, alpha=0.9, lw=0)
        ax.set(xscale="log", yscale="log", xlim=(8e-5, 1.2), ylim=(8e-5, 1.2), xticks=[1e-3, 1e-1], yticks=[1e-4, 1e-2, 1],
               xlabel="p, 999 draws", ylabel="p, 9,999 draws")
        ax.set_xticklabels(["0.001", "0.1"]); ax.set_yticklabels(["0.0001", "0.01", "1"])
        ax.minorticks_off(); clean(ax, grid="both")
        return {"n_components": len(rr), "registered": x, "extended": y}
    return draw


def panel_f(d):
    def draw(ax):
        rows = d.training["published_brain_rows"]
        counts = [d.training["geneformer_source_rows"][0]["cells"], sum(r["cell_number_reserved"] for r in rows)]
        assert counts[1] == 45508
        hbars(ax, ["Geneformer V2", "scFoundation"], counts, [S.BLUE, S.TEAL], xlim=(0, 92000),
              xticks=[0, 30000, 60000], xticklabels=["0", "30k", "60k"], xlabel="Inventory cells", height=0.55,
              annotate=[f"{c:,}" for c in counts], annotate_dx=1500)
        return {"models": ["Geneformer V2", "scFoundation"], "inventory_cells": counts}
    return draw


def panel_g(d):
    levels = d.levels
    exact = [r["exact_cell_ledger_available"] for r in d.training["checkpoint_resolution"]]
    assert not any(exact) and all(row[3] == 0 for row in levels["rows"].values())

    def draw(ax):
        names = list(levels["rows"])
        mat = np.array([levels["rows"][n] for n in names])
        heat(ax, mat, cmap=ListedColormap([S.LIGHT, S.AMBER, S.TEAL]), vmin=0, vmax=2, text=False, row_labels=names,
             col_labels=levels["columns"], col_rotation=90, grid_color="white")
        legend_above(ax, [Patch(facecolor=S.TEAL, label="Verified"), Patch(facecolor=S.AMBER, label="Partial"),
                          Patch(facecolor=S.LIGHT, edgecolor=S.GRAY, lw=0.5, label="Open")],
                     handlelength=0.7, handletextpad=0.25, columnspacing=0.6, bbox_to_anchor=(-0.66, 1.0))
        return {"models": names, "levels": mat.tolist(), "columns": levels["columns"], "scale": levels["scale"]}
    return draw


def panel_h(d):
    def draw(ax):
        rows = d.training["published_brain_rows"]
        selected = set(d.training["selected_sample_ids"])
        vals = [r["cell_number_reserved"] for r in rows]
        colors = [S.TEAL if r["official_sample_id"] in selected else S.ORANGE for r in rows]
        ax.bar(range(len(rows)), vals, color=colors, width=0.72)
        names = [r["official_sample_id"].replace("Sample-", "S") for r in rows]
        ax.set(xlim=(-0.7, len(rows) - 0.3), ylim=(0, 5600), yticks=[0, 2500, 5000], ylabel="Cells")
        ax.set_yticklabels(["0", "2,500", "5,000"])
        ax.set_xticks(range(len(rows))); ax.set_xticklabels(names, rotation=90, fontsize=S.FS_TICK)
        clean(ax, grid="y")
        return {"samples": names, "cells": vals, "selected": [r["official_sample_id"] in selected for r in rows]}
    return draw


def panel_i(d):
    def draw(ax):
        sc = d.training["study_coverage"]
        ns = [sc["training_project_count"], sc["study_project_count"]]
        assert ns == [656, 522] and sc["training_projects_absent_in_study_table"] == 134
        hbars(ax, ["Training", "Study"], ns, [S.TEAL, S.GRAY], xlim=(0, 800), xticks=[0, 300, 600], xlabel="Projects",
              height=0.55, annotate=[str(n) for n in ns], annotate_dx=15)
        return {"training_projects": ns[0], "study_projects": ns[1], "absent_from_study": 134}
    return draw


def build(out_height=130.0):
    d = Data()
    c = Canvas("Figure2", out_height)
    row1, row2, row3 = split(c.rect, "v", [3.1, 4.5, 1.6])
    a_cell, b_cell, e_cell = split(row1, "h", [1.1, 1.0, 1.1])
    c_cell, d_cell, right = split(row2, "h", [1.35, 1.0, 1.1])
    f_cell, g_cell = split(right, "v", [1.0, 1.55])
    h_cell, i_cell = split(row3, "h", [2.2, 1.1])
    c.panel(a_cell, "A", "Incremental alignment", panel_a(d), align_y="top")
    c.panel(b_cell, "B", "Matched co-expression", panel_b(d), align_y="top")
    c.panel(e_cell, "C", "Historical resolution", panel_e(d))
    c.panel(c_cell, "D", "Component probabilities", panel_c(d), align_y="mid")
    c.panel(d_cell, "E", "999 vs 9,999 draws", panel_d(d), align_y="mid")
    c.panel(f_cell, "F", "Source-cell inventory", panel_f(d))
    c.panel(g_cell, "G", "Evidence levels", panel_g(d))
    c.panel(h_cell, "H", "Published brain RNA samples", panel_h(d))
    c.panel(i_cell, "I", "scFoundation catalogues", panel_i(d))
    return c, {}
