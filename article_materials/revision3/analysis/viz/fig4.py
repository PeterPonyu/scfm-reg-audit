"""Figure 4: archived component screen, operational gates and historical readout diagnostics."""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm

from . import style as S
from .common import (FIGDATA, HIST, boolean_heat, clean, fmt, group_labels, heat, load, short)
from .layout import Canvas, split

TEALMAP = LinearSegmentedColormap.from_list("teal", ["#F3F7F7", S.TEAL])
GATE_NAMES = ["Dual full", "Concordance", "Degree sign", "Readout sign", "Above base", "All five"]
METHODS = ["BH", "BY", "Holm", "IUT-BH", "IUT-BY", "IUT-Holm"]


def public(name):
    return load(HIST / "public_results" / f"{name}_v2.public.json")


class Data:
    def __init__(self):
        self.fixed = load(HIST / "analysis/inputs/v2/fixed_panel_audit_v2.json")
        self.rows = [(t, r) for t in ("brain", "pbmc") for r in self.fixed["pooled"][t]["primary_family"]["rows"]]
        self.gate = load(FIGDATA / "gate_external_mapping.json")["rows"]
        self.mult = load(FIGDATA / "multiplicity_robustness.json")["frozen_n999"]["support_counts"]
        self.iut = load(HIST / "analysis/iut_holm_portable.json")["cases"]
        self.attention = public("readout_attention")
        self.ko = public("insilico_ko")
        self.nulls = public("degree_preserving_null")
        self.guard = load(HIST / "analysis/inputs/revision_05613/attention_guard.json")
        assert len(self.rows) == 13 and len(self.gate) == 13
        for (t, r), g in zip(self.rows, self.gate):
            assert (t, r["model_label"]) == (g["tissue"], g["readout"])

    def selection_counts(self):
        out = {}
        for m in METHODS[:5]:
            out[m] = (self.mult[m]["full"], self.mult[m]["non_degree"])
        counts = []
        for spec in ("full", "non_degree"):
            counts.append(sum(c["support_count"] for c in self.iut if c["N"] == 999 and c["spec"] == spec))
        out["IUT-Holm"] = tuple(counts)
        return out


def _row_labels(d):
    return [short(r["model_label"]) for _, r in d.rows]


def _groups():
    return [("Mixed", 0, 7), ("PBMC", 8, 12)]


def panel_a(d):
    def draw(ax):
        m = np.array([[r["mantel"]["p_mc"], r["degree_preserving"]["p_mc"], r["mantel"]["bh_q_family"],
                       r["degree_preserving"]["bh_q_family"]] for _, r in d.rows])
        labels = _row_labels(d)
        heat(ax, m, cmap="Blues", norm=LogNorm(0.001, 1), fmt="{:.3f}", dark=lambda v: v > 0.15, row_labels=labels,
             col_labels=["pM", "pD", "BH\nqM", "BH\nqD"], separators=(8,), col_rotation=0)
        group_labels(ax, _groups(), labels)
        return {"rows": labels, "columns": ["pM", "pD", "BH qM", "BH qD"], "values": m.tolist()}
    return draw


def panel_b(d):
    def draw(ax):
        gp = np.array([g["gate_pass"] for g in d.gate], int)
        full = np.column_stack([gp, gp.all(axis=1)]).astype(float)
        passing = full.sum(axis=0).astype(int)
        omit = [int(np.delete(gp, i, axis=1).all(axis=1).sum()) for i in range(5)]
        assert passing[-1] == 0
        padded = np.vstack([full, np.full((2, 6), np.nan)])
        cmap = ListedColormap([S.LIGHT, S.TEAL]).copy()
        cmap.set_bad("white")
        masked = np.ma.masked_invalid(padded)
        ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
        for k in range(1, 15):
            ax.axhline(k - 0.5, color="white", lw=0.7)
        for k in range(1, 6):
            ax.axvline(k - 0.5, color="white", lw=0.7)
        ax.axhline(7.5, color=S.GRAY, lw=0.8)
        ax.axhline(12.5, color=S.GRAY, lw=0.8)
        for j in range(6):
            ax.text(j, 13, str(passing[j]), ha="center", va="center", fontsize=S.FS_NOTE, fontweight="normal")
            ax.text(j, 14, str(omit[j]) if j < 5 else "–", ha="center", va="center", fontsize=S.FS_NOTE, fontweight="normal")
        labels = _row_labels(d) + ["Pass", "Others pass"]
        ax.set_yticks(range(15)); ax.set_yticklabels(labels)
        ax.set_xticks(range(6)); ax.set_xticklabels(GATE_NAMES, rotation=90)
        ax.tick_params(length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        group_labels(ax, _groups(), labels)
        return {"rows": _row_labels(d), "gates": GATE_NAMES, "pass": full.astype(int).tolist(), "pass_counts": passing.tolist(),
                "omit_one_counts": omit}
    return draw


def panel_c(d):
    def draw(ax):
        x = [r["mantel"]["bh_q_family"] for _, r in d.rows]
        y = [r["degree_preserving"]["bh_q_family"] for _, r in d.rows]
        colors = [S.BLUE if t == "brain" else S.TEAL for t, _ in d.rows]
        ax.scatter(x, y, c=colors, s=17, zorder=3, lw=0)
        ax.axvline(0.05, color=S.RED, ls="--", lw=0.7)
        ax.axhline(0.05, color=S.RED, ls="--", lw=0.7)
        ax.set(xscale="log", yscale="log", xlim=(0.0009, 1.2), ylim=(0.0009, 1.2), xticks=[0.001, 0.05, 1],
               yticks=[0.001, 0.05, 1], xlabel="BH qM", ylabel="BH qD")
        ax.set_xticklabels(["0.001", "0.05", "1"]); ax.set_yticklabels(["0.001", "0.05", "1"])
        ax.minorticks_off(); clean(ax, grid="both")
        return {"qM": x, "qD": y, "context": [t for t, _ in d.rows]}
    return draw


def panel_d(d):
    def draw(ax):
        counts = d.selection_counts()
        mat = np.array([[counts[m][0], counts[m][1]] for m in METHODS], float)
        heat(ax, mat, cmap=TEALMAP, vmin=0, vmax=13, fmt="{:.0f}", dark=lambda v: v > 5, row_labels=METHODS,
             col_labels=["Own", "None"], col_rotation=0)
        return {"methods": METHODS, "own_degree": mat[:, 0].astype(int).tolist(), "no_degree": mat[:, 1].astype(int).tolist(),
                "of_rows": 13}
    return draw


def panel_i(d):
    def draw(ax):
        gr = d.guard["guards"]["without_attention"]["rows"]
        mat = np.array([[r["support"], r["audit_support"]] for r in gr], int)
        assert all(r["agrees"] for r in gr)
        boolean_heat(ax, mat, colors=(S.LIGHT, S.TEAL), row_labels=[short(r["model_label"]) for r in gr],
                     col_labels=["Omit", "Frozen"])
        return {"rows": [short(r["model_label"]) for r in gr], "columns": ["Omit attention", "Frozen screen"], "support": mat.tolist()}
    return draw


def panel_e(d):
    def draw(ax):
        for i, t in enumerate(("brain", "pbmc")):
            for m, off, color, marker, name in (("geneformer_embed", -0.2, S.BLUE, "o", "Embedding"),
                                                ("geneformer_attn", 0.2, S.TEAL, "s", "Attention")):
                r = next(r for tt, r in d.rows if tt == t and r["model_label"] == m)
                ax.scatter(r["observed_partial_rho"], i + off, c=color, marker=marker, s=22, zorder=3, lw=0,
                           label=name if i == 0 else None)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(ylim=(1.7, -0.7), xlim=(-0.005, 0.0065), xticks=[-0.004, 0, 0.004], xlabel="Own-degree ρ")
        ax.set_xticklabels(["−0.004", "0", "0.004"])
        ax.set_yticks([0, 1]); ax.set_yticklabels(["Mixed", "PBMC"])
        ax.legend(loc="center right", handletextpad=0.3, borderaxespad=0.1, handlelength=0.9, labelspacing=0.25)
        clean(ax)
        return {"embedding": [next(r for tt, r in d.rows if tt == t and r["model_label"] == "geneformer_embed")["observed_partial_rho"] for t in ("brain", "pbmc")],
                "attention": [next(r for tt, r in d.rows if tt == t and r["model_label"] == "geneformer_attn")["observed_partial_rho"] for t in ("brain", "pbmc")]}
    return draw


def panel_f(d):
    def draw(ax):
        names = {"brain_geneformer_embed": "GF emb.", "brain_scfoundation": "scFoundation", "brain_uce": "UCE",
                 "brain_coexp_self": "Co-expr."}
        rows = d.nulls
        labels = []
        for i, r in enumerate(rows):
            labels.append(names.get(r["label"], r["label"]))
            ax.errorbar(r["null_mean"], i, xerr=1.96 * r["null_sd"], fmt="o", color=S.GRAY, ms=3.2, capsize=1.8, lw=0.8)
            ax.scatter(r["observed"], i, c=S.TEAL, s=20, zorder=3, lw=0)
        ax.axvline(0, color=S.GRAY, lw=0.5)
        ax.set(ylim=(len(rows) - 0.4, -0.6), xlim=(-0.0035, 0.0068), xticks=[-0.002, 0.002, 0.006],
               xlabel="Observed vs null ρ")
        ax.set_xticklabels(["−0.002", "0.002", "0.006"])
        ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels)
        clean(ax)
        return {"rows": [r["label"] for r in rows], "observed": [r["observed"] for r in rows],
                "null_mean": [r["null_mean"] for r in rows], "null_sd": [r["null_sd"] for r in rows]}
    return draw


def panel_g(d):
    def draw(ax):
        v = d.attention["observed"]
        keys = ["attn_sym_vs_atac", "attn_tf2target_vs_atac", "attn_target2tf_vs_atac", "attn_sym_partial_given_coexp",
                "attn_tf2target_partial_given_coexp"]
        labels = ["Sym.", "TF→target", "Target→TF", "Sym. | cx", "TF→tgt | cx"]
        vals = [v[k] for k in keys]
        ax.barh(range(5), vals, height=0.6, color=[S.BLUE] * 3 + [S.ORANGE] * 2)
        for i, x in enumerate(vals):
            ax.text(x - 0.002, i, fmt(x, 3), ha="right", va="center", fontsize=S.FS_NOTE)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(ylim=(4.55, -0.55), xlim=(-0.092, 0.006), xticks=[-0.04, 0], xlabel="ρ (named design)")
        ax.set_xticklabels(["−0.04", "0"])
        ax.set_yticks(range(5)); ax.set_yticklabels(labels)
        clean(ax)
        return {"designs": labels, "rho": vals}
    return draw


def panel_h(d):
    def draw(ax):
        ko = d.ko
        for off, vals, color, tag in ((-0.2, [ko["ko_vs_atac"], ko["ko_partial_given_coexp"]], S.BLUE, "raw"),
                                      (0.2, [ko["ko_ctrl_vs_atac"], ko["ko_ctrl_partial_given_coexp"]], S.GRAY, "corr.")):
            ax.barh(np.arange(2) + off, vals, height=0.36, color=color)
            for y, x in zip(np.arange(2) + off, vals):
                ax.text(x + 0.003, y, fmt(x, 3), va="center", fontsize=S.FS_NOTE)
                ax.text(0.003, y, tag, va="center", ha="left", fontsize=S.FS_NOTE, color="white")
        ax.set(ylim=(1.6, -0.6), xlim=(0, 0.115), xticks=[0, 0.04, 0.08], xlabel="Historical ρ")
        ax.set_yticks([0, 1]); ax.set_yticklabels(["Marginal", "| co-expr."])
        clean(ax)
        return {"raw_ko": [ko["ko_vs_atac"], ko["ko_partial_given_coexp"]],
                "corrected_ko": [ko["ko_ctrl_vs_atac"], ko["ko_ctrl_partial_given_coexp"]]}
    return draw


def build(out_height=126.0):
    d = Data()
    c = Canvas("Figure4", out_height)
    top, mid = split(c.rect, "v", [86.0, 38.0])
    a_cell, b_cell, right = split(top, "h", [1.0, 1.05, 0.8])
    c_cell, d_cell, i_cell = split(right, "v", [1.25, 1.0, 1.0])
    e_cell, f_cell, g_cell, h_cell = split(mid, "h", [0.85, 1.0, 1.15, 0.9])
    c.panel(a_cell, "A", "Archived component screen", panel_a(d))
    c.panel(b_cell, "B", "Five operational gates", panel_b(d))
    c.panel(c_cell, "C", "BH q concordance", panel_c(d))
    c.panel(d_cell, "D", "Selections by method", panel_d(d))
    c.panel(i_cell, "E", "Attention omission", panel_i(d))
    c.panel(e_cell, "F", "Readout signs", panel_e(d), align_y="r2")
    c.panel(f_cell, "G", "Row-shuffle reference", panel_f(d), align_y="r2")
    c.panel(g_cell, "H", "Attention direction", panel_g(d), align_y="r2")
    c.panel(h_cell, "I", "Knockout control", panel_h(d), align_y="r2")
    return c, {}
