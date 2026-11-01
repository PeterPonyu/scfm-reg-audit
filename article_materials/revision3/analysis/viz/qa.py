"""Geometry and white-space checks computed on the rendered canvas."""
import io
from typing import Dict, List

import numpy as np
from PIL import Image
from matplotlib.text import Text

from . import style as S
from .layout import Canvas

TOL_MM = 0.15
BLOCK_PX = 4


def _text_boxes(canvas: Canvas):
    fig = canvas.fig
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    scale = canvas.w / fig.bbox.width
    boxes = []
    for t in fig.findobj(Text):
        if not (t.get_visible() and t.get_text().strip()):
            continue
        bb = t.get_window_extent(renderer)
        boxes.append((t, (bb.x0 * scale, bb.y0 * scale, bb.x1 * scale, bb.y1 * scale)))
    return boxes


def text_report(canvas: Canvas) -> Dict:
    boxes = _text_boxes(canvas)
    overlaps: List[List[str]] = []
    for i, (t, a) in enumerate(boxes):
        for u, b in boxes[i + 1:]:
            if min(a[2], b[2]) - max(a[0], b[0]) > TOL_MM and min(a[3], b[3]) - max(a[1], b[1]) > TOL_MM:
                overlaps.append([t.get_text(), u.get_text()])
    outside = [t.get_text() for t, a in boxes
               if a[0] < -TOL_MM or a[1] < -TOL_MM or a[2] > canvas.w + TOL_MM or a[3] > canvas.h + TOL_MM]
    out_of_cell = []
    for p in canvas.panels:
        if not p.letter:
            continue
        c = p.cell
        for t in p.ax.findobj(Text):
            if not (t.get_visible() and t.get_text().strip()):
                continue
            bb = t.get_window_extent(canvas.fig.canvas.get_renderer())
            s = canvas.w / canvas.fig.bbox.width
            if bb.x0 * s < c.x0 - 0.3 or bb.x1 * s > c.x1 + 0.3 or bb.y0 * s < c.y0 - 0.3 or bb.y1 * s > c.y1 + 0.3:
                out_of_cell.append([p.letter, t.get_text()])
        for art in (p.title_artist,):
            if art is not None:
                bb = art.get_window_extent(canvas.fig.canvas.get_renderer())
                s = canvas.w / canvas.fig.bbox.width
                if bb.x1 * s > c.x1 + 0.3:
                    out_of_cell.append([p.letter, "title:" + art.get_text()])
    sizes = [t.get_fontsize() for t, _ in boxes]
    return {"text_count": len(boxes), "minimum_font_pt": min(sizes), "overlaps": overlaps,
            "outside_canvas": outside, "outside_cell": out_of_cell}


def _largest_empty_rectangle(white: np.ndarray):
    """Largest all-True axis-aligned rectangle (in cells) by the histogram-stack algorithm."""
    rows, cols = white.shape
    heights = np.zeros(cols, int)
    best = (0, 0, 0, 0, 0)  # area, last row, first col, end col, height
    for r in range(rows):
        heights = np.where(white[r], heights + 1, 0)
        stack: List[tuple] = []  # (start column, height), heights strictly increasing
        for c in range(cols + 1):
            h = int(heights[c]) if c < cols else 0
            start = c
            while stack and stack[-1][1] >= h:
                idx, hh = stack.pop()
                area = hh * (c - idx)
                if area > best[0]:
                    best = (area, r, idx, c, hh)
                start = idx
            stack.append((start, h))
    return best


def blank_bands(white: np.ndarray, px_mm: float):
    """Thickest interior all-white horizontal band and vertical corridor (mm), canvas margins excluded."""
    def longest(flags):
        idx = np.flatnonzero(~flags)
        if idx.size == 0:
            return 0.0
        inner = flags[idx[0]:idx[-1] + 1]
        best = run = 0
        for f in inner:
            run = run + 1 if f else 0
            best = max(best, run)
        return best / px_mm
    return round(longest(white.all(axis=1)), 2), round(longest(white.all(axis=0)), 2)


def render_array(canvas: Canvas, dpi: int = 150) -> np.ndarray:
    buf = io.BytesIO()
    canvas.fig.savefig(buf, format="png", dpi=dpi)
    buf.seek(0)
    return np.asarray(Image.open(buf).convert("RGB"))


def whitespace_report(canvas: Canvas, dpi: int = 150) -> Dict:
    img = render_array(canvas, dpi)
    h, w, _ = img.shape
    px_mm = dpi / 25.4
    white = img.min(axis=2) >= 250
    hb, wb = h // BLOCK_PX, w // BLOCK_PX
    blocks = white[:hb * BLOCK_PX, :wb * BLOCK_PX].reshape(hb, BLOCK_PX, wb, BLOCK_PX).all(axis=(1, 3))
    area, r_end, c0, c1, height = _largest_empty_rectangle(blocks)
    mm_per_block = BLOCK_PX / px_mm
    empty_mm2 = area * mm_per_block ** 2
    canvas_mm2 = canvas.w * canvas.h
    coverage = {}
    for p in canvas.panels:
        if not p.letter:
            continue
        x0, x1 = int(p.cell.x0 * px_mm), int(p.cell.x1 * px_mm)
        y0, y1 = int((canvas.h - p.cell.y1) * px_mm), int((canvas.h - p.cell.y0) * px_mm)
        sub = white[max(y0, 0):y1, max(x0, 0):x1]
        coverage[p.letter] = round(float(1 - sub.mean()), 3)
    band, corridor = blank_bands(white, px_mm)
    x0, x1 = c0 * mm_per_block, c1 * mm_per_block
    y1_top = (r_end - height + 1) * mm_per_block
    y_top, y_bot = y1_top, (r_end + 1) * mm_per_block
    cx, cy = (x0 + x1) / 2, canvas.h - (y_top + y_bot) / 2
    holder = next((p.letter for p in canvas.panels if p.letter and p.cell.x0 <= cx <= p.cell.x1 and p.cell.y0 <= cy <= p.cell.y1), "")
    return {"largest_empty_rectangle_panel": holder, "largest_empty_rectangle_centre_mm": [round(cx, 1), round(cy, 1)],
            "white_fraction": round(float(white.mean()), 3),
            "thickest_blank_row_band_mm": band, "thickest_blank_column_corridor_mm": corridor,
            "largest_empty_rectangle_mm": [round((c1 - c0) * mm_per_block, 1), round(height * mm_per_block, 1)],
            "largest_empty_rectangle_fraction": round(empty_mm2 / canvas_mm2, 4),
            "panel_ink_coverage": coverage, "minimum_panel_ink": min(coverage.values()) if coverage else None}


def gap_report(canvas: Canvas) -> Dict:
    """Smallest and largest clear distance between neighbouring inner axes rectangles (mm)."""
    rects = [(p.letter, p.ax.get_position()) for p in canvas.panels if p.letter]
    rects = [(k, (b.x0 * canvas.w, b.y0 * canvas.h, b.x1 * canvas.w, b.y1 * canvas.h)) for k, b in rects]
    overlaps = []
    for i, (ka, a) in enumerate(rects):
        for kb, b in rects[i + 1:]:
            if min(a[2], b[2]) - max(a[0], b[0]) > 0.05 and min(a[3], b[3]) - max(a[1], b[1]) > 0.05:
                overlaps.append([ka, kb])
    return {"inner_axes_overlaps": overlaps}


def full_report(canvas: Canvas) -> Dict:
    report = {"figure": canvas.name, "size_mm": [round(canvas.w, 2), round(canvas.h, 2)]}
    report.update(text_report(canvas))
    report.update(whitespace_report(canvas))
    report.update(gap_report(canvas))
    report["panels"] = [p.letter for p in canvas.panels if p.letter]
    report["inner_axes_mm"] = canvas.inner_rects()
    # Directly verify the requested typography in the rendered artist tree.
    labels = {id(p.letter_artist) for p in canvas.panels if p.letter_artist is not None}
    report["non_label_bold_text"] = [t.get_text() for t, _ in _text_boxes(canvas)
                                      if str(t.get_fontweight()) in {"bold", "heavy", "semibold", "700", "800", "900"}
                                      and id(t) not in labels]
    report["title_style_errors"] = []
    for p in canvas.panels:
        if p.title_artist is None:
            continue
        title = p.title_artist
        error = abs(title.get_position()[0] * canvas.w - (p.cell.x0 + p.cell.x1) / 2)
        if title.get_ha() != "center" or str(title.get_fontweight()) != "normal" or error > 1e-6:
            report["title_style_errors"].append(p.letter)
    report["title_count"] = sum(p.title_artist is not None for p in canvas.panels)
    report["passed"] = (not report["overlaps"] and not report["outside_canvas"] and not report["outside_cell"]
                        and not report["inner_axes_overlaps"] and not report["non_label_bold_text"]
                        and not report["title_style_errors"] and report["minimum_font_pt"] >= S.MIN_FONT_PT - 1e-6)
    return report
