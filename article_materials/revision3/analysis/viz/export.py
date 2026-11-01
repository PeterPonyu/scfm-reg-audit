"""Write PDF, editable SVG with embedded Arial, PNG preview, 600-dpi JPEG and the data record."""
import base64
import io
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from fontTools import subset
from fontTools.ttLib import TTFont

from . import style as S
from .layout import Canvas

SVG_NS = "http://www.w3.org/2000/svg"


def _embed_svg_fonts(path: Path):
    """Keep SVG text editable and embed subset WOFF Arial so browsers never substitute a family."""
    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree = ET.parse(path)
    root = tree.getroot()
    text = "".join("".join(n.itertext()) for n in root.iter() if n.tag in (f"{{{SVG_NS}}}text", f"{{{SVG_NS}}}tspan"))
    rules = []
    for weight, key in (("normal", "regular"), ("bold", "bold")):
        font = TTFont(str(S.FONT_PATHS[key]), recalcTimestamp=False)
        sub = subset.Subsetter(options=subset.Options(recalc_timestamp=False))
        sub.populate(text=text)
        sub.subset(font)
        font.flavor = "woff"
        buf = io.BytesIO()
        font.save(buf)
        data = base64.b64encode(buf.getvalue()).decode()
        rules.append("@font-face{font-family:Arial;font-weight:%s;font-style:normal;"
                     "src:url(data:font/woff;base64,%s) format(\"woff\");}" % (weight, data))
    defs = root.find(f"{{{SVG_NS}}}defs")
    if defs is None:
        defs = ET.SubElement(root, f"{{{SVG_NS}}}defs")
    ET.SubElement(defs, f"{{{SVG_NS}}}style", {"type": "text/css"}).text = "\n".join(rules)
    tree.write(path, encoding="utf-8", xml_declaration=True)
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


def export(canvas: Canvas, out_dir: Path, stem: str, data: dict, report: dict, jpeg_dpi: int = 600):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = canvas.fig
    fig.savefig(out_dir / f"{stem}.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(out_dir / f"{stem}.svg", metadata={"Date": None})
    _embed_svg_fonts(out_dir / f"{stem}.svg")
    fig.savefig(out_dir / f"{stem}.png", dpi=300)
    fig.savefig(out_dir / f"{stem}.jpg", dpi=jpeg_dpi, pil_kwargs={"quality": 95})
    record = {"figure": stem, "display_only": True, "panels": canvas.records(), "layout": report,
              "typography": {"family": "Arial", "minimum_pt": S.MIN_FONT_PT},
              "encodings": {"no_degree": S.ORANGE, "own_degree": S.TEAL, "reference_degree": S.BLUE}}
    record.update(data)
    (out_dir / f"{stem}.data.json").write_text(json.dumps(record, indent=1, default=_json_default) + "\n")


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))
