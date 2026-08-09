#!/usr/bin/env python3
"""Rebuild the PeerJ submission package from the canonical manuscript.

Generates:
- source/: self-contained editable article source (TikZ fragments included)
- internal/figure_build/: standalone wrappers for each figure
- flat_upload/: flat manuscript with pre-rendered FigureN.pdf files, tables,
  references, and the compiled manuscript PDF
- SHA256SUMS.txt for every flat_upload file

Typography/layout contract:
- Font: 11pt newtxtext/newtxmath (Times-compatible, matches PeerJ NimbusRomNo9L)
- Engine: pdftex (tikzDevice compatibility)
- Design width: 6.8in -> ~6.5in in manuscript
- Unified with make_figs.R via figure_typography.py
"""
import hashlib
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
PKG = PAPER / "submission_peerj"
TABLE_FRAGMENTS = [
    "table1_primary_fixed_panel",
    "table2_cross_tissue_observed",
    "table3_pertype_ranges",
    "table4_protocol_pass",
    "table5_related_work",
]

# Import typography contract (figure map, validation)
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_typography import FIGURE_MAP, validate_figure_map, TypographyContract

FIG_NAMES = list(FIGURE_MAP)
TYPOGRAPHY = TypographyContract()

WRAPPER = """\\documentclass[tikz,border=0pt]{standalone}
\\usepackage{amsmath}
\\usepackage{newtxtext}
\\usepackage{newtxmath}
\\begin{document}
\\input{FIGFILE}
\\end{document}
"""


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _check_log(cwd, pdf_name):
    log_path = Path(cwd) / (Path(pdf_name).stem + ".log")
    if log_path.exists():
        log_text = log_path.read_text(errors="replace")
        for bad in ("undefined citation", "undefined reference",
                    "Label(s) may have changed", "There were undefined citations"):
            if bad.lower() in log_text.lower():
                raise RuntimeError(f"LaTeX log contains '{bad}' in {log_path}")


def build_pdf_with_bib(cwd, tex_name="manuscript.tex"):
    """Build a PDF that has a bibliography using an explicit, deterministic
    pdflatex -> bibtex -> pdflatex -> pdflatex sequence.

    latexmk's automatic bibtex handling is unreliable on a clean tree in some
    TeX Live builds: the first pdflatex pass returns a nonzero exit on the
    (expected) undefined-citation warnings and latexmk halts before running
    bibtex, leaving citations unresolved. Driving the passes explicitly avoids
    that ordering hazard. pdflatex may still return nonzero on its final pass
    due to benign warnings, so success is judged by the log gate below and the
    presence of the PDF, not by the exit code."""
    stem = Path(tex_name).stem
    for stale in (".aux", ".bbl", ".blg"):
        p = Path(cwd) / f"{stem}{stale}"
        if p.exists():
            p.unlink()
    pdflatex = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_name]
    subprocess.run(pdflatex, cwd=cwd, capture_output=True, text=True)
    subprocess.run(["bibtex", stem], cwd=cwd, capture_output=True, text=True)
    subprocess.run(pdflatex, cwd=cwd, capture_output=True, text=True)
    last = subprocess.run(pdflatex, cwd=cwd, capture_output=True, text=True)
    # A third pdflatex can still leave "Label(s) may have changed"; one more
    # pass is required before the log gate, otherwise packaging aborts on a
    # clean tree after figure regeneration.
    log_path = Path(cwd) / f"{stem}.log"
    if log_path.exists() and "Label(s) may have changed" in log_path.read_text(errors="replace"):
        last = subprocess.run(pdflatex, cwd=cwd, capture_output=True, text=True)
    pdf_path = Path(cwd) / f"{stem}.pdf"
    if not pdf_path.exists():
        raise RuntimeError(
            f"build failed: no {stem}.pdf produced in {cwd}\n"
            f"{last.stdout[-2000:]}\n{last.stderr[-2000:]}")
    _check_log(cwd, f"{stem}.pdf")
    return last


def run(cmd, cwd, pdf_name=None):
    """Run a build command; latexmk may exit non-zero on a clean first pass when the
    bbl does not exist yet, while still producing a valid PDF. Retry once in that case.
    Afterwards, reject any log carrying undefined citations/references."""
    attempts = 2 if pdf_name else 1
    last = None
    for _ in range(attempts):
        last = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        if last.returncode == 0:
            break
        if pdf_name and (Path(cwd) / pdf_name).exists():
            continue
        break
    if last.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\n{last.stdout[-2000:]}\n{last.stderr[-2000:]}")
    if pdf_name:
        _check_log(cwd, pdf_name)
    return last


def _copy_fig_assets(name, dest_dir):
    """Copy a figure fragment's companion raster assets (e.g. rasterized colorbars
    emitted by tikzDevice as {name}_ras*.png) next to the .tex so \\pgfimage/
    \\includegraphics can resolve them at compile time."""
    for asset in sorted((PAPER / "figs").glob(f"{name}_ras*.png")):
        shutil.copy2(asset, Path(dest_dir) / asset.name)


def build_source():
    source = PKG / "source"
    if source.exists():
        shutil.rmtree(source)
    (source / "figs").mkdir(parents=True)
    shutil.copy2(PAPER / "manuscript.tex", source / "manuscript.tex")
    shutil.copy2(PAPER / "references.bib", source / "references.bib")
    shutil.copy2(PAPER / "wlpeerj.cls", source / "wlpeerj.cls")
    for name, _ in FIG_NAMES:
        shutil.copy2(PAPER / "figs" / f"{name}.tex", source / "figs" / f"{name}.tex")
        _copy_fig_assets(name, source / "figs")
    for name in TABLE_FRAGMENTS:
        shutil.copy2(PAPER / "figs" / f"{name}.tex", source / "figs" / f"{name}.tex")
    build_pdf_with_bib(source, "manuscript.tex")
    return source


def build_figures():
    wrappers = PKG / "internal" / "figure_build"
    if wrappers.exists():
        shutil.rmtree(wrappers)
    wrappers.mkdir(parents=True)
    flat = PKG / "flat_upload"
    flat.mkdir(parents=True, exist_ok=True)
    for name, figure in FIG_NAMES:
        (wrappers / f"{figure}.tex").write_text(WRAPPER.replace("FIGFILE", f"{name}.tex"))
        shutil.copy2(PAPER / "figs" / f"{name}.tex", wrappers / f"{name}.tex")
        _copy_fig_assets(name, wrappers)
        run(["latexmk", "-pdf", "-interaction=nonstopmode", f"{figure}.tex"], wrappers,
            pdf_name=f"{figure}.pdf")
        shutil.copy2(wrappers / f"{figure}.pdf", flat / f"{figure}.pdf")
    return flat


def build_flat(flat):
    for stale in flat.glob("manuscript.*"):
        if stale.suffix not in {".tex", ".pdf"}:
            stale.unlink()
    text = (PAPER / "manuscript.tex").read_text()
    text = text.replace("\\usepackage{tikz}\n", "\\usepackage{xcolor}\n")
    text = text.replace("\\graphicspath{{figs/}}\n", "")
    for name, figure in FIG_NAMES:
        text = text.replace(
            f"\\fitfig{{\\input{{figs/{name}.tex}}}}",
            f"\\includegraphics[width=\\linewidth]{{{figure}.pdf}}",
        )
        text = re.sub(
            rf"\\resizebox\{{([^}}]+)\}}\{{!\}}\{{\\input\{{figs/{re.escape(name)}\.tex\}}\}}",
            rf"\\includegraphics[width=\1]{{{figure}.pdf}}",
            text,
        )
    text = text.replace("figs/table1_primary_fixed_panel.tex", "table1_primary_fixed_panel.tex")
    text = text.replace("figs/table2_cross_tissue_observed.tex", "table2_cross_tissue_observed.tex")
    text = text.replace("figs/table3_pertype_ranges.tex", "table3_pertype_ranges.tex")
    text = text.replace("figs/table4_protocol_pass.tex", "table4_protocol_pass.tex")
    text = text.replace("figs/table5_related_work.tex", "table5_related_work.tex")
    (flat / "manuscript.tex").write_text(text)
    shutil.copy2(PAPER / "references.bib", flat / "references.bib")
    shutil.copy2(PAPER / "wlpeerj.cls", flat / "wlpeerj.cls")
    for name in TABLE_FRAGMENTS:
        shutil.copy2(PAPER / "figs" / f"{name}.tex", flat / f"{name}.tex")
    build_pdf_with_bib(flat, "manuscript.tex")
    for stale in flat.iterdir():
        if stale.is_file() and stale.suffix in {".aux", ".blg", ".fdb_latexmk",
                                                 ".fls", ".log", ".out"}:
            stale.unlink()
    return flat


def write_checksums(flat):
    lines = []
    for path in sorted(flat.iterdir()):
        if path.is_file() and path.suffix in {".pdf", ".tex", ".bib", ".cls"}:
            lines.append(f"{sha256_file(path)}  {path.name}")
    (PKG / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
    return lines


def build_upload_zip(flat, checksum_names):
    """Build a clean, deterministic upload archive containing ONLY the files
    listed in SHA256SUMS.txt (plus the checksum file itself). Building from the
    manifest rather than ``zip -r flat_upload`` guarantees that transient
    tooling state (.omc/, .pytest_cache/, editor scratch) and LaTeX
    intermediates can never leak into the archive that a human uploads."""
    import zipfile
    zip_path = PKG / "upload.zip"
    if zip_path.exists():
        zip_path.unlink()
    members = sorted(checksum_names) + ["SHA256SUMS.txt"]
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in members:
            src = (flat / name) if name != "SHA256SUMS.txt" else (PKG / name)
            if not src.exists():
                raise RuntimeError(f"upload manifest lists missing file: {name}")
            zf.write(src, arcname=name)
    return zip_path


def regenerate_figures():
    """Regenerate panel data and canonical TeX before packaging."""
    run(["python3", "make_panel_data.py"], PAPER)
    run(["Rscript", "make_figs.R"], PAPER)


def scrub_generated_files():
    """Remove caches and LaTeX intermediates from submission directories."""
    cache_names = {".mypy_cache", ".pytest_cache", ".ruff_cache", ".omc"}
    build_dirs = (PKG / "source", PKG / "flat_upload", PKG / "internal" / "figure_build")
    for directory in build_dirs:
        if not directory.exists():
            continue
        for path in list(directory.rglob("*")):
            if path.is_dir() and path.name in cache_names:
                shutil.rmtree(path)
    for path in (PKG / "source").glob("manuscript.*"):
        if path.suffix not in {".tex", ".bbl", ".pdf"}:
            path.unlink()
    for path in (PKG / "internal" / "figure_build").glob("*"):
        if path.is_file() and path.suffix in {".aux", ".blg", ".fdb_latexmk", ".fls", ".log", ".out"}:
            path.unlink()


def refresh_canonical_pdf():
    """Rebuild paper/manuscript.pdf so the canonical human-review PDF never
    drifts from the freshly regenerated figures and the packaged mirror.
    Without this, editing a figure and repackaging leaves paper/manuscript.pdf
    stale, which previously shipped an outdated figure caption/title into the
    canonical review target."""
    build_pdf_with_bib(PAPER, "manuscript.tex")
    for stale in PAPER.glob("manuscript.*"):
        if stale.suffix not in {".tex", ".bbl", ".pdf"}:
            stale.unlink()
    rplots = PAPER / "Rplots.pdf"
    if rplots.exists():
        rplots.unlink()


def main():
    regenerate_figures()
    validate_figure_map(str(PAPER / "manuscript.tex"))
    source = build_source()
    flat = build_figures()
    build_flat(flat)
    refresh_canonical_pdf()
    scrub_generated_files()
    if list(PKG.rglob(".omc")):
        raise RuntimeError(".omc session state present in submission package")
    lines = write_checksums(flat)
    checksum_names = [ln.split("  ", 1)[1] for ln in lines]
    zip_path = build_upload_zip(flat, checksum_names)
    print(f"source rebuilt: {source}")
    print(f"flat upload rebuilt: {flat} ({len(lines)} files)")
    print(f"clean upload archive: {zip_path} ({len(checksum_names) + 1} members)")
    print(f"canonical PDF refreshed: {PAPER / 'manuscript.pdf'}")


if __name__ == "__main__":
    main()
