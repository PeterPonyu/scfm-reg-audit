"""Figure 1 (R3): design, source correction, registered admission and common support."""
import sys

import numpy as np

from . import style as S
from .common import ROOT, clean, short
from .layout import Canvas, fixed, split
from . import fig1 as old
from .schematics import proximal, setup, brain, cells, route, matrix, graph, text

sys.path.insert(0, str(ROOT / "analysis"))
import r3lib as L  # noqa: E402

CTX = {"brain_selected_panel-brain": "New / brain", "brain_selected_panel-pbmc": "New / PBMC",
       "historical_panel_corrected_source-brain": "Old / brain"}


def workflow(ax):
    setup(ax)
    brain(ax, 11, 83, 11, S.GRAY)
    brain(ax, 44, 83, 11, S.BLUE)
    cells(ax, 78, 83, 11)
    text(ax, 16.5, 72, "Aorta RNA +\nbrain ATAC", color=S.GRAY)
    text(ax, 50, 72, "Genuine brain\nfull pool, halves", color=S.BLUE)
    text(ax, 83.5, 72, "BMMC\nfive donors", color=S.TEAL)
    for x in (16.5, 50, 83.5):
        route(ax, [(x, 63), (x, 59), (50, 59)], arrow=False)
    route(ax, [(50, 59), (50, 55)])
    matrix(ax, 8, 44, 11, 9)
    text(ax, 60, 48, "Panel → common support")
    route(ax, [(50, 42), (50, 38), (25, 38), (25, 34)])
    route(ax, [(50, 38), (77, 38), (77, 34)])
    graph(ax, 8, 23, 14, 10, S.BLUE)
    graph(ax, 61, 23, 14, 10, S.TEAL)
    text(ax, 38, 28, "RNA graph", color=S.BLUE)
    text(ax, 87, 28, "ATAC\nproxy", color=S.TEAL)
    route(ax, [(25, 21), (25, 17), (50, 17)], arrow=False)
    route(ax, [(77, 21), (77, 17), (50, 17), (50, 13)])
    text(ax, 50, 8, "Degree designs → M/D → IUT-BY / Holm")
    return {"schematic": True, "sources": ["aorta RNA / brain ATAC audit record", "genuine brain full pool and disjoint halves",
                                           "BMMC five evaluation donors"]}


class Data:
    def __init__(self):
        self.old = old.Data()
        self.support = {}
        for job in L.JOBS:
            ctx = L.load_context(job)
            six = L.common_support(ctx)[4]
            gl = L.common_support(ctx, pair_readouts=L.GENE_LEVEL)[4]
            self.support[job] = {"six": six, "gene_level": gl}


def panel_g(d):
    def draw(ax):
        metrics = [("genes", "Genes"), ("tfs", "TFs"), ("pairs", "Pairs, six readouts"), ("pairs_gl", "Pairs, five gene-level")]
        colors = [S.BLUE, S.PURPLE, S.TEAL, S.SAND]
        rec, y, ticks, names = [], 0.0, [], []
        for job in L.JOBS:
            s, g = d.support[job]["six"], d.support[job]["gene_level"]
            vals = {"genes": s["genes_S"] / s["genes_panel"], "tfs": s["tfs_S"] / s["tfs_panel"],
                    "pairs": s["fraction_registered_pairs_retained"], "pairs_gl": g["fraction_registered_pairs_retained"]}
            for k, ((key, lab), col) in enumerate(zip(metrics, colors)):
                yy = y + k * 0.8
                ax.barh(yy, 100 * vals[key], height=0.66, color=col, lw=0)
                ax.text(100 * vals[key] + 1.2, yy, f"{100 * vals[key]:.1f}", va="center", fontsize=S.FS_NOTE)
            ticks.append(y + 1.2)
            names.append(CTX[job])
            rec.append({"context": CTX[job], **{k: float(v) for k, v in vals.items()}})
            y += 4.0
        ax.set(xlim=(0, 118), xticks=[0, 50, 100], ylim=(y - 0.6, -0.8), xlabel="Retained in common support (%)")
        ax.set_yticks(ticks)
        ax.set_yticklabels(names)
        ax.tick_params(axis="y", length=0)
        clean(ax)
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=c, label=lab) for (k, lab), c in zip(metrics, colors)], loc="upper center",
                  bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, handlelength=1.0, fontsize=S.FS_NOTE,
                  columnspacing=0.9, handletextpad=0.4)
        return {"rows": rec, "measured_slots": {"registered": 16, "common_support": 72}}
    return draw


def build(out_height=150.0):
    d = Data()
    c = Canvas("Figure1", out_height)
    schematic_row, results = split(c.rect, "v", [fixed(48.0), 1.0])
    workflow_cell, proximal_cell = split(schematic_row, "h", [1.0, 1.0])
    c.panel(workflow_cell, "A", "Source-to-audit workflow", workflow)
    c.panel(proximal_cell, "B", "Proximal proxy construction", proximal)
    top, bot = split(results, "v", [1.0, 1.55])
    c_cell, d_cell, e_cell = split(top, "h", [1.0, 0.95, 1.1])
    f_cell, g_cell = split(bot, "h", [1.15, 1.0])
    c.panel(c_cell, "C", "Source fingerprints", old.panel_a(d.old), align_y="r1")
    c.panel(d_cell, "D", "Panel overlap", old.panel_c(d.old), align_y="r1")
    c.panel(e_cell, "E", "Selection quota", old.panel_b(d.old), align_y="r1")
    c.panel(f_cell, "F", "Registered gene / pair coverage (%)", old.panel_g(d.old, legend=True), align_y="r2")
    c.panel(g_cell, "G", "Common support", panel_g(d), align_y="r2")
    return c, {"support": d.support}
