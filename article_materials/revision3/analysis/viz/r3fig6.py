"""Figure 6 (R3): disjoint discovery and evaluation samples, and documented training exposure."""
import json

import numpy as np
from matplotlib.patches import Rectangle

from . import style as S
from .common import ROOT, clean, legend_above, pt, short
from .layout import Canvas, split
from . import fig2 as old

RES = ROOT / "evidence" / "R3_RESULTS.json"
ORDER = ["geneformer_embed", "geneformer_attention", "random_init_geneformer", "scgpt_encoder", "scfoundation_encoder", "uce_encoder"]
DEG = {"none": (S.ORANGE, "^"), "own": (S.TEAL, "s")}


class Data:
    def __init__(self):
        r = json.loads(RES.read_text())
        self.split = r["split"]
        assert self.split["streams_found"] == 8
        self.rows = self.split["rows"]
        self.cs = r["common_support"]["six"]["rows"]
        self.old = old.Data()
        # Weight-identity evidence for scFoundation from this revision's checksum and tensor checks.
        scf = json.loads((ROOT / "evidence" / "scfoundation_verification.json").read_text())
        rows = self.old.levels["rows"]
        key = next(k for k in rows if k.lower().startswith("scfoundation"))
        col = next(i for i, c in enumerate(self.old.levels["columns"]) if c.lower().startswith("weight"))
        # A match to the GenBio-hosted bundle alone does not establish developer-release identity.
        level = (2 if scf.get("developer_origin_verified") and scf.get("local_identical_to_official_submodel")
                 else (1 if scf.get("local_equals_mirror_file") else 0))
        rows[key][col] = level
        self.scf_weight_level = level

    def pick(self, fold, m, deg, null="restricted"):
        return next(x for x in self.rows if (x["fold"], x["readout"], x["degree_mode"], x["null_mode"]) == (fold, m, deg, null))


def assignment_panel(d):
    def draw(ax):
        a = d.split["assignment"]
        rec = {}
        for row, half in enumerate(("half_A", "half_B")):
            for k, s in enumerate(sorted(a[half], key=lambda v: (a["diagnosis"][v], a["batch"][v], v))):
                col = S.RED if a["diagnosis"][s] == "AD" else S.BLUE
                ax.add_patch(Rectangle((k, row * 1.3), 0.92, 0.9, facecolor=col, edgecolor="none"))
                ax.text(k + 0.46, row * 1.3 + 0.45, s.replace("Sample-", ""), ha="center", va="center", color="white",
                        fontsize=S.FS_NOTE)
                ax.text(k + 0.46, row * 1.3 - 0.18, f"b{a['batch'][s]}", ha="center", va="center", fontsize=S.FS_NOTE - 0.5)
            ax.text(-0.25, row * 1.3 + 0.45, "Half A" if half == "half_A" else "Half B", ha="right", va="center",
                    fontsize=S.FS_NOTE)
            rec[half] = a[half]
        ax.set(xlim=(-2.2, 9.1), ylim=(2.4, -0.45))
        ax.axis("off")
        legend_above(ax, [pt("s", S.RED, "AD"), pt("s", S.BLUE, "Control")])
        return {"halves": rec}
    return draw


def overlap_panel(d):
    def draw(ax):
        prep = d.split["prepared"]
        rows = []
        for f in ("1", "2"):
            p = prep[f]
            rows.append((f"Fold {f}", p["overlap_with_registered_brain_panel"], p["panel_tfs"]))
        y = np.arange(len(rows))
        ax.barh(y, [r[1] for r in rows], height=0.55, color=S.TEAL)
        ax.barh(y, [1200 - r[1] for r in rows], left=[r[1] for r in rows], height=0.55, color=S.LIGHT, edgecolor=S.GRAY, lw=0.4)
        for yy, r in zip(y, rows):
            ax.text(r[1] / 2, yy, f"{r[1]:,} shared", ha="center", va="center", color="white", fontsize=S.FS_NOTE)
            ax.text(1210, yy, f"{r[2]} TFs", va="center", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, 1400), xticks=[0, 600, 1200], ylim=(1.6, -0.6), xlabel="Discovery-panel genes")
        ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
        clean(ax)
        return {"rows": [{"fold": r[0], "shared_with_full_pool_panel": r[1], "tfs": r[2]} for r in rows]}
    return draw


def effect_panel(d):
    def draw(ax):
        rec = []
        for k, m in enumerate(ORDER):
            for deg, off in (("none", -0.18), ("own", 0.18)):
                col, mk = DEG[deg]
                v1, v2 = d.pick(1, m, deg)["observed"], d.pick(2, m, deg)["observed"]
                full = next(x for x in d.cs if (x["job"], x["readout"], x["degree_mode"], x["null_mode"]) ==
                            ("brain_selected_panel-brain", m, deg, "restricted"))["observed"]
                ax.plot([v1, v2], [k + off, k + off], color=col, lw=0.7, zorder=2)
                ax.scatter([v1], [k + off], marker=mk, color=col, s=13, lw=0, zorder=3)
                ax.scatter([v2], [k + off], marker=mk, facecolor="white", edgecolor=col, s=13, lw=0.8, zorder=3)
                ax.scatter([full], [k + off], marker="|", color="black", s=26, lw=0.9, zorder=4)
                rec.append({"readout": m, "degree": deg, "fold1": v1, "fold2": v2, "full_pool_common_support": full})
        ax.axvline(0, color=S.GRAY, lw=0.6)
        ax.set(xlim=(-0.012, 0.037), xticks=[-0.01, 0, 0.01, 0.02, 0.03], ylim=(5.6, -0.7),
               xlabel="Incremental partial ρ, evaluation half")
        ax.set_xticklabels(["−0.01", "0", "0.01", "0.02", "0.03"])
        ax.set_yticks(range(6)); ax.set_yticklabels([short(m) for m in ORDER])
        clean(ax)
        legend_above(ax, [pt("o", S.GRAY, "Fold 1"), pt("o", S.GRAY, "Fold 2", hollow=True), pt("|", "black", "Full pool")])
        return {"rows": rec}
    return draw


def pvalue_panel(d):
    def draw(ax):
        rec = []
        for k, m in enumerate(ORDER):
            for fold, off in ((1, -0.2), (2, 0.2)):
                for deg in ("none", "own"):
                    r = d.pick(fold, m, deg)
                    col, mk = DEG[deg]
                    ax.scatter(-np.log10(r["pIUT"]), k + off, marker=mk, s=12, lw=0.7, zorder=3,
                               facecolor=col if fold == 1 else "white", edgecolor=col)
                    rec.append({"fold": fold, "readout": m, "degree": deg, "pIUT": r["pIUT"], "q_BY": r["q_BY"], "p_Holm": r["p_Holm"]})
        m48 = 48
        ax.axvline(-np.log10(0.05 / m48), color=S.RED, ls="--", lw=0.7)
        ax.axvline(-np.log10(0.05), color=S.GRAY, ls=":", lw=0.7)
        ax.set(xlim=(0, 4.3), xticks=[0, 1, 2, 3, 4], ylim=(5.6, -0.7), xlabel="−log10 joint p, restricted")
        ax.set_yticks(range(6)); ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        clean(ax)
        return {"rows": rec, "family_slots": m48}
    return draw


def build(out_height=112.0):
    d = Data()
    c = Canvas("Figure6", out_height)
    top, mid = split(c.rect, "v", [0.72, 1.35])
    a_cell, b_cell = split(top, "h", [1.6, 1.0])
    c_cell, d_cell, e_cell = split(mid, "h", [1.25, 0.8, 1.15])
    c.panel(a_cell, "A", "Stratified sample halves", assignment_panel(d))
    c.panel(b_cell, "B", "Panel stability", overlap_panel(d))
    c.panel(c_cell, "C", "Held-out effects", effect_panel(d), align_y="m")
    c.panel(d_cell, "D", "Joint p", pvalue_panel(d), align_y="m")
    c.panel(e_cell, "E", "Training-exposure evidence", old.panel_g(d.old), align_y="m")
    return c, {"split_family": d.split.get("family_48")}
