"""Millimetre layout engine for dense, non-grid multi-panel figures.

A figure is partitioned by nested ``split`` calls into outer cells. Every cell receives
one panel. A panel is drawn once, then *fitted*: its axes are moved and resized so that
the tight bounding box of tick labels, axis labels and annotations (the outer box) fills
the cell minus the header line. Gaps between panels are therefore exactly the explicit
``GAP_*`` constants, and panels in one alignment group share their inner extent.
"""
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.text import Text

from . import style as S

HEADER_MM = 3.7          # letter + short title line above each panel
GAP_X = 1.7              # explicit horizontal gap between outer cells
GAP_Y = 1.5              # explicit vertical gap between outer cells
EDGE_MM = 0.6            # canvas edge padding on every side
LETTER_ADVANCE_MM = 4.3  # minimum clearance reserved for the panel letter


@dataclass(frozen=True)
class Rect:
    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def w(self):
        return self.x1 - self.x0

    @property
    def h(self):
        return self.y1 - self.y0

    def inset(self, l=0.0, b=0.0, r=0.0, t=0.0):
        return Rect(self.x0 + l, self.y0 + b, self.x1 - r, self.y1 - t)

    def as_list(self):
        return [round(v, 3) for v in (self.x0, self.y0, self.x1, self.y1)]


Fixed = Tuple[str, float]
Size = Union[float, Fixed]


def fixed(mm: float) -> Fixed:
    return ("mm", float(mm))


def split(rect: Rect, axis: str, sizes: Sequence[Size], gap: Optional[float] = None) -> List[Rect]:
    """Split ``rect``. axis 'h' runs left to right, 'v' runs top to bottom.

    Each size is a relative weight (float) or ``fixed(mm)``.
    """
    if gap is None:
        gap = GAP_X if axis == "h" else GAP_Y
    total = rect.w if axis == "h" else rect.h
    free = total - gap * (len(sizes) - 1)
    fixed_total = sum(s[1] for s in sizes if isinstance(s, tuple))
    weight_total = sum(s for s in sizes if not isinstance(s, tuple))
    lengths = [s[1] if isinstance(s, tuple) else (free - fixed_total) * s / weight_total for s in sizes]
    out, pos = [], (rect.x0 if axis == "h" else rect.y1)
    for length in lengths:
        if axis == "h":
            out.append(Rect(pos, rect.y0, pos + length, rect.y1))
            pos += length + gap
        else:
            out.append(Rect(rect.x0, pos - length, rect.x1, pos))
            pos -= length + gap
    return out


@dataclass
class _Panel:
    letter: str
    title: str
    cell: Rect
    ax: plt.Axes
    record: dict
    header: bool
    align_left: Optional[str]
    align_y: Optional[str]
    pad: Tuple[float, float, float, float]   # extra clearance (l, b, r, t) in mm
    nudge: Tuple[float, float, float, float]  # manual (dx0, dy0, dx1, dy1) on the final inner rect
    fixed_inner: Optional[Rect] = None
    letter_artist: Optional[Text] = None
    title_artist: Optional[Text] = None


class Canvas:
    def __init__(self, name: str, height_mm: float, width_mm: float = S.TEXT_WIDTH_MM):
        self.name, self.w, self.h = name, width_mm, height_mm
        self.fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), dpi=100)
        self.panels: List[_Panel] = []
        self.extra_texts: List[Text] = []

    # ---- geometry helpers -------------------------------------------------
    @property
    def rect(self) -> Rect:
        return Rect(EDGE_MM, EDGE_MM, self.w - EDGE_MM, self.h - EDGE_MM)

    def _frac(self, r: Rect):
        return [r.x0 / self.w, r.y0 / self.h, r.w / self.w, r.h / self.h]

    def mm_to_fig(self, x, y):
        return x / self.w, y / self.h

    # ---- panels -----------------------------------------------------------
    def panel(self, cell: Rect, letter: str, title: str, draw: Callable[[plt.Axes], Optional[dict]], *,
              header: bool = True, align_left: Optional[str] = None, align_y: Optional[str] = None,
              pad=(0.0, 0.0, 0.0, 0.0), nudge=(0.0, 0.0, 0.0, 0.0)) -> plt.Axes:
        top = HEADER_MM if header else 0.0
        seed = cell.inset(l=10, b=7, r=1, t=top + 1)
        if seed.w < 4 or seed.h < 4:  # tiny strips: provisional box only; fitting resizes it
            seed = Rect(cell.x0 + 1, cell.y0 + 1, max(cell.x1 - 1, cell.x0 + 6), max(cell.y1 - top - 1, cell.y0 + 6))
        ax = self.fig.add_axes(self._frac(seed))
        record = draw(ax) or {}
        p = _Panel(letter, title, cell, ax, record, header, align_left, align_y, tuple(pad), tuple(nudge))
        if header:
            lx, ly = self.mm_to_fig(cell.x0 + 0.2, cell.y1 - 0.35)
            p.letter_artist = self.fig.text(lx, ly, letter, fontsize=S.FS_LETTER, fontweight="bold",
                                            ha="left", va="top", gid=f"panel-{letter}-label")
            if title:
                tx, ty = self.mm_to_fig((cell.x0 + cell.x1) / 2, cell.y1 - 0.75)
                p.title_artist = self.fig.text(tx, ty, title, fontsize=S.FS_TITLE, fontweight="normal",
                                               ha="center", va="top", gid=f"panel-{letter}-title")
        self.panels.append(p)
        return ax

    def note(self, cell: Rect, text: str, *, size=S.FS_NOTE, color="black", ha="left", va="center",
             weight="normal", linespacing=1.15, bbox=None) -> Text:
        x = {"left": cell.x0, "center": (cell.x0 + cell.x1) / 2, "right": cell.x1}[ha]
        y = {"top": cell.y1, "center": (cell.y0 + cell.y1) / 2, "bottom": cell.y0}[va]
        fx, fy = self.mm_to_fig(x, y)
        t = self.fig.text(fx, fy, text, fontsize=size, color=color, ha=ha, va=va, fontweight=weight,
                          linespacing=linespacing, bbox=bbox)
        self.extra_texts.append(t)
        return t

    def key(self, cell: Rect, handles, ncol: int, *, columnspacing=1.2, handlelength=1.1) -> plt.Axes:
        """Legend strip occupying a whole cell (no header)."""
        ax = self.fig.add_axes(self._frac(cell))
        ax.axis("off")
        ax.legend(handles=handles, loc="center", ncol=ncol, frameon=False, handlelength=handlelength,
                  handletextpad=0.45, columnspacing=columnspacing, borderpad=0, borderaxespad=0,
                  labelspacing=0.3, fontsize=S.FS_NOTE)
        self.panels.append(_Panel("", "", cell, ax, {}, False, None, None, (0, 0, 0, 0), (0, 0, 0, 0), fixed_inner=cell))
        return ax

    # ---- fitting ----------------------------------------------------------
    def _decorations(self, p: _Panel, renderer):
        """Return (left, bottom, right, top) clearance of the outer box around the inner axes, in mm."""
        fig_bb = self.fig.bbox
        tb = p.ax.get_tightbbox(renderer)
        ab = p.ax.get_window_extent(renderer)
        px_x, px_y = self.w / fig_bb.width, self.h / fig_bb.height
        return ((ab.x0 - tb.x0) * px_x, (ab.y0 - tb.y0) * px_y, (tb.x1 - ab.x1) * px_x, (tb.y1 - ab.y1) * px_y)

    def fit(self, iterations: int = 3):
        targets = [p for p in self.panels if p.fixed_inner is None]
        for _ in range(iterations):
            self.fig.canvas.draw()
            renderer = self.fig.canvas.get_renderer()
            deco = {id(p): self._decorations(p, renderer) for p in targets}
            for p in targets:
                l, b, r, t = deco[id(p)]
                if p.align_y:
                    peers = [q for q in targets if q.align_y == p.align_y]
                    b, t = max(deco[id(q)][1] for q in peers), max(deco[id(q)][3] for q in peers)
                if p.align_left:
                    peers = [q for q in targets if q.align_left == p.align_left]
                    l = max(deco[id(q)][0] for q in peers)
                top = HEADER_MM if p.header else 0.0
                avail = p.cell.inset(t=top)
                inner = Rect(avail.x0 + l + p.pad[0], avail.y0 + b + p.pad[1], avail.x1 - r - p.pad[2],
                             avail.y1 - t - p.pad[3])
                nx0, ny0, nx1, ny1 = p.nudge
                inner = Rect(inner.x0 + nx0, inner.y0 + ny0, inner.x1 + nx1, inner.y1 + ny1)
                p.ax.set_position(self._frac(inner))
        self.fig.canvas.draw()

    # ---- output -----------------------------------------------------------
    def inner_rects(self) -> Dict[str, List[float]]:
        out = {}
        for p in self.panels:
            bb = p.ax.get_position()
            out[p.letter or "key"] = [round(bb.x0 * self.w, 3), round(bb.y0 * self.h, 3),
                                      round(bb.x1 * self.w, 3), round(bb.y1 * self.h, 3)]
        return out

    def records(self) -> dict:
        return {p.letter: p.record for p in self.panels if p.letter}
