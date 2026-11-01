"""Figure 2 (R3): calibration of unrestricted and restricted component tests under synthetic stress."""
import json

import numpy as np

from . import style as S
from .common import ROOT, clean, legend_above, pt
from .layout import Canvas, split

CAL = ROOT / "evidence" / "calibration"
CONTEXTS = [("brain", "Mixed-source (historical)"), ("pbmc", "PBMC (historical panel)"), ("genuine_brain", "Genuine brain")]
SCEN = [("iid", "Independent"), ("shock_stratifier_scaled", "Shocks: strata"), ("shock_detection_scaled", "Shocks: detection")]
COMP = [("M", S.BLUE, "o", "Gene label M"), ("D", S.ORANGE, "s", "Target shuffle D"), ("joint", S.TEAL, "D", "Joint (IUT)")]
NULLS = [("unrestricted", "U"), ("restricted", "R")]


class Data:
    def __init__(self):
        self.rep = {c: json.loads((CAL / c / "calibration.json").read_text()) for c, _ in CONTEXTS}
        self.arr = {(c, n): dict(np.load(CAL / c / f"{n}.npz", allow_pickle=False)) for c, _ in CONTEXTS for n, _ in NULLS}
        for c, _ in CONTEXTS:
            assert self.rep[c]["status"] == "complete"


def row(d, ctx, null, scen):
    return next(r for r in d.rep[ctx][null] if r["scenario"] == scen)


def rate_panel(d, ctx, labels):
    def draw(ax):
        rec, ticks, names, y = [], [], [], 0.0
        ax.axvspan(0.025, 0.075, color=S.LIGHT, lw=0, zorder=0)
        for s, sname in SCEN:
            for n, nshort in NULLS:
                r = row(d, ctx, n, s)
                for j, (k, col, mk, _) in enumerate(COMP):
                    v = r[k]
                    lo, hi = v["wilson95"]
                    yy = y + (j - 1) * 0.24
                    ax.errorbar(v["rate"], yy, xerr=[[v["rate"] - lo], [hi - v["rate"]]], fmt=mk, color=col, ms=3.0,
                                mfc=col if n == "restricted" else "white", capsize=1.4, lw=0.75, zorder=3)
                    rec.append({"scenario": s, "null": n, "component": k, "count": v["count"], "rate": v["rate"],
                                "wilson95": v["wilson95"], "nominal": v["nominal_bradley_and_wilson"]})
                jt = r["joint"]
                ax.text(0.468, y, f"{jt['count']}", ha="right", va="center", fontsize=S.FS_NOTE, color=S.TEAL)
                ticks.append(y)
                names.append(f"{sname} {nshort}")
                y += 1.0
            y += 0.45
        ax.axvline(0.05, color=S.RED, ls="--", lw=0.75, zorder=2)
        ax.text(0.468, -0.95, "Joint/200", ha="right", va="center", fontsize=S.FS_NOTE, color=S.TEAL)
        ax.set(xlim=(0, 0.48), ylim=(y - 0.25, -1.3), xticks=[0, 0.1, 0.2, 0.3, 0.4], xlabel="Rejection rate at p ≤ 0.05")
        ax.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4"])
        ax.set_yticks(ticks)
        ax.set_yticklabels(names if labels else [])
        if not labels:
            ax.tick_params(axis="y", length=0)
        clean(ax)
        return {"context": ctx, "rows": rec}
    return draw


def ecdf_panel(d, scen):
    def draw(ax):
        grid = np.linspace(0, 0.2, 401)
        ax.plot([0, 0.2], [0, 0.2], color=S.GRAY, lw=0.7, ls=":")
        rec = []
        colors = {"brain": S.PURPLE, "pbmc": S.BLUE, "genuine_brain": S.TEAL}
        for c, cname in CONTEXTS:
            for n, nshort in NULLS:
                a = d.arr[(c, n)]
                sel = a["labels"] == scen
                p = np.maximum(a["pM"][sel], a["pD"][sel])
                f = (p[None, :] <= grid[:, None]).mean(1)
                ax.step(grid, f, where="post", color=colors[c], lw=0.95 if n == "restricted" else 0.8,
                        ls="-" if n == "restricted" else (0, (2.2, 1.4)))
                rec.append({"context": c, "null": n, "fraction_le_0.05": float((p <= 0.05).mean()), "n": int(sel.sum())})
        ax.axvline(0.05, color=S.RED, ls="--", lw=0.6)
        ax.set(xlim=(0, 0.2), ylim=(0, 0.42), xticks=[0, 0.05, 0.1, 0.15, 0.2], yticks=[0, 0.1, 0.2, 0.3, 0.4],
               xlabel="Joint p threshold", ylabel="Fraction rejected")
        ax.set_xticklabels(["0", "0.05", "0.1", "0.15", "0.2"])
        clean(ax, grid="y")
        return {"scenario": scen, "rows": rec}
    return draw


def build(out_height=104.0):
    d = Data()
    c = Canvas("Figure2", out_height)
    key, top, bottom = split(c.rect, "v", [0.09, 1.0, 0.78])
    a_cell, b_cell, c_cell = split(top, "h", [1.42, 1.0, 1.0])
    d_cell, e_cell, f_cell = split(bottom, "h", [1.0, 1.0, 0.9])
    handles = [pt(m, col, n) for _, col, m, n in COMP] + [pt("o", S.GRAY, "U hollow / R filled")]
    c.key(key, handles, ncol=4)
    c.panel(a_cell, "A", "Mixed-source (historical)", rate_panel(d, "brain", True), align_y="top")
    c.panel(b_cell, "B", "Old-panel PBMC (historical)", rate_panel(d, "pbmc", False), align_y="top")
    c.panel(c_cell, "C", "New panel / genuine brain", rate_panel(d, "genuine_brain", False), align_y="top")
    c.panel(d_cell, "D", "Joint p, shocks on strata", ecdf_panel(d, "shock_stratifier_scaled"), align_y="bot")
    c.panel(e_cell, "E", "Joint p, shocks on detection", ecdf_panel(d, "shock_detection_scaled"), align_y="bot")

    def legend(ax):
        ax.axis("off")
        from matplotlib.lines import Line2D
        hs = [Line2D([], [], color=S.PURPLE, lw=1.0, label="Mixed-source"), Line2D([], [], color=S.BLUE, lw=1.0, label="Old-panel PBMC"),
              Line2D([], [], color=S.TEAL, lw=1.0, label="Genuine brain"),
              Line2D([], [], color="black", lw=0.95, label="Restricted"),
              Line2D([], [], color="black", lw=0.8, ls=(0, (2.2, 1.4)), label="Unrestricted"),
              Line2D([], [], color=S.GRAY, lw=0.7, ls=":", label="Uniform")]
        ax.legend(handles=hs, loc="center left", frameon=False, fontsize=S.FS_NOTE, handlelength=1.8, labelspacing=0.45)
        return {}
    c.panel(f_cell, "F", "Key for D and E", legend, header=True)
    return c, {"source": "evidence/calibration/*/calibration.json and *.npz"}
