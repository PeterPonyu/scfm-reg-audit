"""Shared data access and drawing helpers for the panel modules."""
import json
from pathlib import Path

import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

from . import style as S

ROOT = Path(__file__).resolve().parents[2]
HIST = ROOT / "historical_evidence"
FIGDATA = HIST / "figures" / "data"
PORTABLE = HIST / "figures" / "corrected" / "portable" / "plotted_data"
RESOLUTION = ROOT / "analysis" / "resolution" / "results.json"
if not RESOLUTION.exists():  # layout of the public deposit
    RESOLUTION = ROOT.parent / "resolution_9999" / "results" / "results.json"

READOUTS = ["geneformer_embed", "geneformer_attention", "random_init_geneformer", "scgpt_encoder",
            "scfoundation_encoder", "uce_encoder"]
ADMITTED = [READOUTS[i] for i in (0, 2, 3, 5)]
JOBS = ["brain_selected_panel-brain", "brain_selected_panel-pbmc", "historical_panel_corrected_source-brain"]
SHORT = {
    "geneformer_embed": "GF emb.", "geneformer_embedding": "GF emb.", "geneformer_attention": "GF attn.",
    "geneformer_attn": "GF attn.", "random_init_geneformer": "Random GF", "random_init_floor": "Random GF",
    "random_floor": "Random GF", "scgpt_encoder": "scGPT", "scGPT_encoder": "scGPT",
    "scfoundation_encoder": "scFoundation", "scFoundation_encoder": "scFoundation", "uce_encoder": "UCE",
    "UCE_encoder": "UCE", "geneformer_ko_raw": "GF KO raw", "geneformer_ko_posctrl": "GF KO corr.",
    "co_expression": "Co-expr.",
}
DEGREE = [("none", S.ORANGE, "^", "No degree"), ("own", S.TEAL, "s", "Own degree")]


def load(path):
    return json.loads(Path(path).read_text())


def payload(name):
    return load(PORTABLE / f"{name}.json")["payload"]


def short(m):
    return SHORT.get(m, m)


def corrected_rows():
    """Registered post hoc brain rows that were measured (16 slots)."""
    return [r for r in load(RESOLUTION)["rows"] if r["measured"]]


def clean(ax, grid="x", ticks=True):
    ax.tick_params(length=2.2 if ticks else 0, width=0.6, pad=1.8)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis=grid, color=S.GRID, lw=0.5)
    return ax


def bar_labels(ax, ys, xs, labels, dx, *, color="black", ha="left", weight="normal"):
    for y, x, t in zip(ys, xs, labels):
        ax.text(x + dx, y, t, va="center", ha=ha, fontsize=S.FS_NOTE, color=color, fontweight=weight)


def heat(ax, matrix, *, cmap, norm=None, vmin=None, vmax=None, fmt="{:.3f}", dark=None, row_labels=None,
         col_labels=None, col_rotation=0, labels_top=False, separators=(), text=True, fontsize=S.FS_NOTE,
         grid_color="white"):
    """Annotated heatmap with white cell gutters. ``dark(value)`` decides white versus black text."""
    matrix = np.asarray(matrix, float)
    ax.imshow(matrix, cmap=cmap, norm=norm, vmin=vmin if norm is None else None, vmax=vmax if norm is None else None,
              aspect="auto", interpolation="nearest")
    rows, cols = matrix.shape
    if text:
        for i in range(rows):
            for j in range(cols):
                v = matrix[i, j]
                if np.isnan(v):
                    continue
                color = "white" if (dark(v) if dark else False) else "black"
                ax.text(j, i, fmt(v) if callable(fmt) else fmt.format(v), ha="center", va="center", fontsize=fontsize, color=color)
    for k in range(1, rows):
        ax.axhline(k - 0.5, color=grid_color, lw=0.7)
    for k in range(1, cols):
        ax.axvline(k - 0.5, color=grid_color, lw=0.7)
    for k in separators:
        ax.axhline(k - 0.5, color=S.GRAY, lw=0.8)
    ax.set_yticks(range(rows))
    ax.set_yticklabels(row_labels if row_labels is not None else [], fontsize=S.FS_TICK)
    ax.set_xticks(range(cols))
    ax.set_xticklabels(col_labels if col_labels is not None else [], fontsize=S.FS_TICK, rotation=col_rotation,
                       ha="center" if col_rotation in (0, 90) else "right",
                       rotation_mode="anchor" if col_rotation not in (0, 90) else "default")
    if labels_top:
        ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return ax


def boolean_heat(ax, matrix, *, colors=(S.LIGHT, S.TEAL), **kwargs):
    kwargs.setdefault("text", False)
    return heat(ax, matrix, cmap=ListedColormap(list(colors)), vmin=0, vmax=len(colors) - 1, **kwargs)


def hbars(ax, labels, values, colors, *, xlim, xticks, xlabel, xticklabels=None, height=0.62, annotate=None,
          annotate_dx=0.0, inside=False, show_labels=True, log=False, ylim=None, left=None):
    """Horizontal bar panel with optional value annotation (outside, or inside in white)."""
    n = len(labels)
    ys = np.arange(n)
    ax.barh(ys, values, height=height, color=colors, left=left)
    if annotate is not None:
        for y, v, t in zip(ys, values, annotate):
            if inside:
                ax.text(v / 2 if not log else np.sqrt(v), y, t, ha="center", va="center", color="white", fontsize=S.FS_NOTE)
            else:
                ax.text(v * (1 + annotate_dx) if log else v + annotate_dx, y, t, ha="left", va="center", fontsize=S.FS_NOTE)
    ax.set(xlim=xlim, ylim=ylim or (n - 0.45, -0.55), xticks=xticks, xlabel=xlabel)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels if show_labels else [])
    if not show_labels:
        ax.tick_params(axis="y", length=0)
    clean(ax)
    return ax


def fmt(v, nd=4, sign=False):
    """Number text with a true minus sign."""
    t = f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"
    return t.replace("-", "\u2212")


def text_width_pt(texts, size=None, weight="normal"):
    """Widest rendered width (points) of the given strings in the active font."""
    import matplotlib as mpl
    from matplotlib.font_manager import FontProperties
    from matplotlib.backends.backend_agg import RendererAgg
    size = size or S.FS_TICK
    renderer = RendererAgg(100, 100, 100)
    prop = FontProperties(family=mpl.rcParams["font.sans-serif"][0], size=size, weight=weight)
    return max(renderer.get_text_width_height_descent(t, prop, ismath=False)[0] for t in texts) * 72 / 100


def group_labels(ax, groups, tick_labels, *, gap_pt=4.0):
    """Rotated group labels left of the y tick labels. groups: (label, first_row, last_row)."""
    offset = -(text_width_pt(tick_labels) + 2.2 + 1.8 + gap_pt + S.FS_NOTE / 2)
    for label, a, b in groups:
        ax.annotate(label, xy=(0, (a + b) / 2), xycoords=("axes fraction", "data"), xytext=(offset, 0),
                    textcoords="offset points", rotation=90, ha="center", va="center", fontsize=S.FS_NOTE,
                    fontweight="normal", annotation_clip=False)
        ax.plot([0, 0], [a - 0.38, b + 0.38], transform=ax.get_yaxis_transform(), color=S.GRAY, lw=0.8,
                solid_capstyle="butt", clip_on=False)


def fmt_rho(v, nd=3):
    """Signed-neutral rounding: values that round to zero print as 0.000 (no negative zero)."""
    r = round(float(v), nd)
    return fmt(0.0 if r == 0 else r, nd)


def legend_above(ax, handles, ncol=None, **overrides):
    """Compact legend directly above the axes; its extent is part of the panel's outer box."""
    options = dict(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=ncol or len(handles), frameon=False, handletextpad=0.3,
                   columnspacing=0.9, borderaxespad=0.15, borderpad=0.0, handlelength=0.9, labelspacing=0.2,
                   fontsize=S.FS_NOTE)
    options.update(overrides)
    return ax.legend(handles=handles, **options)


def pt(marker, color, label, hollow=False, ms=4.4):
    return Line2D([], [], marker=marker, ls="none", color=color, mfc="white" if hollow else color, ms=ms, label=label)
