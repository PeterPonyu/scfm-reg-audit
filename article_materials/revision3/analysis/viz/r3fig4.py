"""Figure 4 (R3): ten degree designs in the three source-matched contexts."""
import json

import numpy as np
from matplotlib.colors import TwoSlopeNorm

from . import style as S
from .common import ROOT, clean, heat, legend_above, short
from .layout import Canvas, split

RES = ROOT / "evidence" / "R3_RESULTS.json"
JOBS = ["brain_selected_panel-brain", "brain_selected_panel-pbmc", "historical_panel_corrected_source-brain"]
CTX = {"brain_selected_panel-brain": "New/brain", "brain_selected_panel-pbmc": "New/PBMC",
       "historical_panel_corrected_source-brain": "Old/brain"}
ORDER = ["geneformer_embed", "geneformer_attention", "random_init_geneformer", "scgpt_encoder", "scfoundation_encoder", "uce_encoder"]
DESIGNS = ["none", "own", "outdeg_only", "indeg_only", "log_degree", "rank_degree", "spline_degree", "strength",
           "external_fibro", "external_other"]
DLAB = {"none": "None", "own": "Own linear", "outdeg_only": "Out only", "indeg_only": "In only", "log_degree": "Log",
        "rank_degree": "Rank", "spline_degree": "Spline", "strength": "Weighted", "external_fibro": "Fibroblast",
        "external_other": "Other tissue"}
# Historical IUT-BY selections out of 13 per design (Revision 2, Table 5).
HISTORICAL = {"none": 0, "own": 7, "outdeg_only": 1, "indeg_only": 0, "log_degree": 1, "rank_degree": 1,
              "spline_degree": 6, "strength": 0, "external_fibro": 0, "external_other": 0}


class Data:
    def __init__(self):
        g = json.loads(RES.read_text())["grid"]
        assert g["designs_found"] == 30
        self.rows = g["rows"]
        self.per = g["per_context_design"]
        self.by_design = g["selection_by_design"]
        self.pooled = g["pooled_180"]

    def get(self, job, m, des):
        return next(r for r in self.rows if (r["job"], r["readout"], r["design"]) == (job, m, des))


def heat_panel(d):
    def draw(ax):
        mat = 1000 * np.array([[d.get(j, m, des)["observed"] for des in DESIGNS] for j in JOBS for m in ORDER])
        labels = [f"{CTX[j]} {short(m)}" for j in JOBS for m in ORDER]
        heat(ax, mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-30, vcenter=0, vmax=30), fmt="{:.0f}",
             dark=lambda v: abs(v) > 18, row_labels=labels, col_labels=[DLAB[x] for x in DESIGNS], col_rotation=90)
        for t in ax.texts:
            v = t.get_text()
            t.set_text("0" if v in ("-0", "−0") else v.replace("-", "−"))
        for k in (6, 12):
            ax.axhline(k - 0.5, color="black", lw=0.8)
        rec = []
        for i, j in enumerate(JOBS):
            for k, m in enumerate(ORDER):
                for c, des in enumerate(DESIGNS):
                    r = d.get(j, m, des)
                    rec.append({"context": CTX[j], "readout": m, "design": des, "rho": r["observed"], "pIUT": r["pIUT"]})
        return {"rows": rec}
    return draw


def count_panel(d):
    def draw(ax):
        y = np.arange(len(DESIGNS))
        hist = [100 * HISTORICAL[x] / 13 for x in DESIGNS]
        new_by = [100 * d.by_design[x]["BY"] / 18 for x in DESIGNS]
        ax.barh(y - 0.19, hist, height=0.34, color=S.GRAY, label="Historical, IUT-BY of 13")
        ax.barh(y + 0.19, new_by, height=0.34, color=S.TEAL, label="Source-matched, IUT-BY of 18")
        for yy, h, b, des in zip(y, hist, new_by, DESIGNS):
            if HISTORICAL[des]:
                ax.text(h + 1, yy - 0.19, str(HISTORICAL[des]), va="center", fontsize=S.FS_NOTE)
            by, ho = d.by_design[des]["BY"], d.by_design[des]["Holm"]
            if by or ho:
                ax.text(b + 1, yy + 0.19, f"{by}/{ho}", va="center", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, 70), xticks=[0, 20, 40, 60], ylim=(len(DESIGNS) - 0.5, -0.6), xlabel="Rows selected (%)")
        ax.set_yticks(y); ax.set_yticklabels([DLAB[x] for x in DESIGNS])
        clean(ax)
        ax.legend(loc="lower right", frameon=False, fontsize=S.FS_NOTE, handlelength=1.0, borderaxespad=0.2)
        return {"historical": HISTORICAL, "source_matched": d.by_design}
    return draw


def width_panel(d):
    def draw(ax):
        rec = []
        for i, j in enumerate(JOBS):
            ratios = [d.get(j, m, "none")["nullM_sd"] / d.get(j, m, "own")["nullM_sd"] for m in ORDER]
            ax.scatter(ratios, np.arange(6) + (i - 1) * 0.22, s=12, color=[S.BLUE, S.TEAL, S.PURPLE][i], lw=0, zorder=3,
                       label=CTX[j])
            rec += [{"context": CTX[j], "readout": m, "ratio": r} for m, r in zip(ORDER, ratios)]
        ax.axvline(1, color=S.GRAY, ls=":", lw=0.7)
        ax.set(xscale="log", xlim=(0.9, 25), ylim=(5.6, -0.7), xlabel="Null SD ratio, no / own degree (M)")
        ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(["1", "2", "5", "10", "20"])
        ax.xaxis.set_minor_locator(__import__("matplotlib").ticker.NullLocator())
        ax.set_yticks(range(6)); ax.set_yticklabels([short(m) for m in ORDER])
        clean(ax)
        legend_above(ax, [__import__("viz.common", fromlist=["pt"]).pt("o", c, CTX[j]) for c, j in zip([S.BLUE, S.TEAL, S.PURPLE], JOBS)])
        return {"rows": rec}
    return draw


def build(out_height=120.0):
    d = Data()
    c = Canvas("Figure4", out_height)
    left, right = split(c.rect, "h", [1.25, 1.0])
    rt, rb = split(right, "v", [1.25, 1.0])
    c.panel(left, "A", "Partial ρ under ten designs", heat_panel(d))
    c.panel(rt, "B", "Selections by design", count_panel(d))
    c.panel(rb, "C", "Gene-label null width", width_panel(d))
    return c, {"pooled_180": d.pooled, "per_context_design": d.per}
