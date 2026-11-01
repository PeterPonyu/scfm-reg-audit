"""Typography, palette and size constants shared by every figure.

All text is Arial (regular and bold only; no italic, no fallback family). Sizes are
final print sizes in points because figures are drawn at the exact manuscript text
width and embedded without scaling.
"""
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "screg-mpl"))
import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

TEXT_WIDTH_MM = 159.2          # A4 text width with 1-inch margins
PT_MM = 25.4 / 72.0

# Final-size typography (points).
FS_TICK = 7.5
FS_LABEL = 7.5
FS_NOTE = 7.5
FS_TITLE = 8.0
FS_LETTER = 10.0
MIN_FONT_PT = 7.0

# Semantic palette. Degree encodings are identical in every figure.
BLUE = "#286B8C"
TEAL = "#237E75"
ORANGE = "#B46D2B"
GRAY = "#8B9299"
LIGHT = "#EDF0F2"
RED = "#A54949"
AMBER = "#F1E4D7"
PURPLE = "#6D75A5"
SAND = "#B4A76C"
GRID = "#E3E7EA"
RULE = "#B8C2C8"
CELLTYPES = [BLUE, TEAL, ORANGE, RED, GRAY, PURPLE, SAND]

FONT_DIR_ENV = "ARIAL_FONT_DIR"


def _configured_font_dir():
    """Font folder from analysis/paths.local.json ("font_dir") if set, else the user font folder."""
    cfg = Path(__file__).resolve().parents[1] / "paths.local.json"
    value = json.loads(cfg.read_text()).get("font_dir") if cfg.exists() else None
    return Path(value).expanduser() if value else Path.home() / ".local/share/fonts"


DEFAULT_FONT_DIR = _configured_font_dir()
FONT_FILES = {"regular": "Arial.ttf", "bold": "Arial Bold.ttf"}


def font_paths():
    base = Path(os.environ.get(FONT_DIR_ENV, DEFAULT_FONT_DIR))
    paths = {k: base / v for k, v in FONT_FILES.items()}
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Arial files not found; set {FONT_DIR_ENV}: {missing}")
    return paths


def apply_typography():
    """Register Arial and set rcParams. Returns the font file paths used."""
    paths = font_paths()
    for p in paths.values():
        font_manager.fontManager.addfont(str(p))
    family = "Arial"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [family],
        "mathtext.fontset": "custom", "mathtext.rm": family, "mathtext.it": family,
        "mathtext.bf": f"{family}:bold", "mathtext.sf": family, "mathtext.default": "regular",
        "mathtext.fallback": None,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.weight": "normal", "axes.titleweight": "normal", "axes.labelweight": "normal",
        "font.size": FS_NOTE, "axes.labelsize": FS_LABEL, "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK, "legend.fontsize": FS_NOTE, "axes.titlesize": FS_TITLE,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2, "xtick.major.pad": 1.8,
        "ytick.major.pad": 1.8, "axes.labelpad": 1.6, "axes.spines.top": False,
        "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.facecolor": "white", "lines.linewidth": 0.9, "legend.frameon": False,
        "svg.hashsalt": "array05613-r2-final", "hatch.linewidth": 0.5,
    })
    return paths


FONT_PATHS = apply_typography()
plt.ioff()
