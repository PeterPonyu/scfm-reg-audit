#!/usr/bin/env python3
"""Checks for the Hugo Pages assemble output. Fail closed."""
from __future__ import annotations

import http.server
import re
import shutil
import socketserver
import sys
import tempfile
import threading
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
PUBLIC = Path(sys.argv[1]) if len(sys.argv) > 1 else SITE / "public"
LAYOUTS = SITE / "layouts"

BANNED = re.compile(
    r"peerj|frontiers|\bgenetics\b|\bmanuscript\b|\bjournal\b|\babstract\b|"
    r"\bsubmission\b|\bvenue\b|\bpublisher\b|\bpreprint\b|\barticle\b|"
    r"computer\s+science|\bphd\b|ph\.d|\bpapers?\b",
    flags=re.I,
)

# Printed order 1-13 then A1-A3. fig10_coverage_qc is printed 13, not 1.
FIG_STEMS = [
    "fig_study_design",
    "fig1_truth_construct",
    "fig2_cross_tissue_decomp",
    "fig3_primary_audit",
    "fig4_usability_check",
    "fig5_null_diagnostics",
    "fig6_spec_sensitivity",
    "fig7_pertype_descriptive",
    "fig8_injection_ladder",
    "fig9_tf_probe",
    "fig11_third_tissue_transfer",
    "fig12_protocol_pass_matrix",
    "fig10_coverage_qc",
    "fig_ext1_construct_mantel",
    "fig_ext2_baselines_collectri",
    "fig_ext3_honesty_policy",
]


def fail(msg: str) -> None:
    print(f"test: FAIL {msg}", file=sys.stderr)
    raise SystemExit(1)


def ok(msg: str) -> None:
    print(f"test: ok {msg}")


def main() -> None:
    if not PUBLIC.is_dir():
        fail(f"missing publish dir {PUBLIC}")

    routes = [
        "index.html",
        "figures/index.html",
        "reproducibility/index.html",
        "404.html",
        ".nojekyll",
    ]
    for rel in routes:
        path = PUBLIC / rel
        if not path.is_file():
            fail(f"missing route {rel}")
        ok(rel)

    for rel in ("peerj/index.html", "frontiers/index.html"):
        if (PUBLIC / rel).is_file():
            fail(f"retired route still published: {rel}")
    ok("retired venue routes absent")

    hub = (PUBLIC / "index.html").read_text(encoding="utf-8")
    if "0/13" not in hub:
        fail("hub missing protocol-pass 0/13")
    ok("hub has 0/13")

    if "446" not in hub or "1,200" not in hub:
        fail("hub missing frozen 446 × 1,200 panel")
    if "Monte Carlo" not in hub or "degree-preserving" not in hub:
        fail("hub missing dual-null wording")
    if "Support" not in hub:
        fail("hub missing Support")
    ok("hub has panel and dual-null")

    if re.search(r"Ph\.D|PhD", hub, flags=re.I):
        fail("hub contains Ph.D / PhD chrome")
    if "Army Medical" in hub or "word-count" in hub.lower():
        fail("hub contains degree / word-count chrome")
    if BANNED.search(hub):
        fail(f"hub contains banned word: {BANNED.search(hub).group(0)!r}")
    if "products/" in hub or ".pdf" in hub.lower():
        fail("hub advertises PDFs")
    ok("hub identity chrome checks")

    for html_path in PUBLIC.rglob("*.html"):
        text = html_path.read_text(encoding="utf-8")
        rel = html_path.relative_to(PUBLIC).as_posix()
        if re.search(r"Ph\.D|PhD", text):
            fail(f"{rel} contains Ph.D / PhD")
        if "Army Medical" in text:
            fail(f"{rel} contains Army Medical")
        if "fonts.googleapis.com" in text or "fonts.gstatic.com" in text:
            fail(f"{rel} loads Google Fonts CDN")
        hit = BANNED.search(text)
        if hit:
            fail(f"{rel} contains banned word: {hit.group(0)!r}")
    ok("no banned words / Ph.D / Army Medical / webfont CDN in HTML")

    for layout in LAYOUTS.rglob("*.html"):
        text = layout.read_text(encoding="utf-8")
        if "/scfm-reg-audit/" in text:
            fail(f"{layout.relative_to(SITE)} hardcodes /scfm-reg-audit/ (use relURL)")
    if 'href="/scfm-reg-audit/' in (LAYOUTS / "index.html").read_text(encoding="utf-8"):
        fail("hub layout hardcodes project-pages prefix")
    ok("layouts do not hardcode /scfm-reg-audit/")

    if "/scfm-reg-audit/figures/" not in hub:
        fail("built hub missing relURL prefix /scfm-reg-audit/figures/")
    if "/scfm-reg-audit/reproducibility/" not in hub:
        fail("built hub missing relURL prefix /scfm-reg-audit/reproducibility/")
    if 'href="/figures/' in hub or 'href="/reproducibility/' in hub:
        fail("built hub has root-absolute science route (missing project base)")
    if re.search(r"href=/>", hub) or 'href="/"' in hub:
        fail("built hub has root-absolute href=/ (missing project base)")
    if 'href=/scfm-reg-audit/' not in hub.replace('href="/scfm-reg-audit/', "href=/scfm-reg-audit/"):
        fail("built hub missing project-pages href prefix")
    if "/peerj/" in hub or "/frontiers/" in hub:
        fail("built hub still links retired routes")
    ok("built hub uses project-pages prefix")

    pdfs = list(PUBLIC.rglob("*.pdf"))
    if pdfs:
        fail(f"expected no PDFs, found {sorted(p.name for p in pdfs)}")
    for name in (f"Figure{i}.pdf" for i in range(1, 13)):
        if (PUBLIC / name).exists() or list(PUBLIC.rglob(name)):
            fail(f"denylist PDF present: {name}")
    ok("publish dir has no PDFs")

    figures = (PUBLIC / "figures/index.html").read_text(encoding="utf-8")
    if "fig10_coverage_qc" not in figures:
        fail("figures page missing coverage QC stem")
    if ">13<" not in figures and ">13</" not in figures:
        if not re.search(r">\s*13\s*<", figures):
            fail("figures page missing printed 13")
    if "fig_ext1_construct_mantel" not in figures:
        fail("figures page missing A1 stem")
    if "paper/" in figures or "submission_" in figures:
        fail("figures page leaks repo path words")
    ok("figures catalog is stem + printed order")

    # The page must show the graphs, not just name them: one visible <img>
    # per stem, src under /scfm-reg-audit/figures/, file present in public/.
    srcs = re.findall(r'<img\b[^>]*?\bsrc="?([^"\s>]+)"?', figures)
    want_srcs = {f"/scfm-reg-audit/figures/{stem}.png" for stem in FIG_STEMS}
    missing_src = sorted(want_srcs - set(srcs))
    if missing_src:
        fail(f"figures page missing <img> for {missing_src}")
    if len([s for s in srcs if s in want_srcs]) != len(FIG_STEMS):
        fail(f"figures page should embed {len(FIG_STEMS)} previews, found {len(srcs)}")
    for stem in FIG_STEMS:
        preview = PUBLIC / "figures" / f"{stem}.png"
        if not preview.is_file():
            fail(f"missing published preview figures/{stem}.png")
        if preview.stat().st_size < 20_000:
            fail(f"preview figures/{stem}.png suspiciously small")
    ok("figures page embeds 16 visible graph previews, all files published")

    smoke(PUBLIC)
    print("test: all checks passed")


def smoke(public: Path) -> None:
    tmp = Path(tempfile.mkdtemp(prefix="scfm-pages-smoke-"))
    mounted = tmp / "scfm-reg-audit"
    shutil.copytree(public, mounted)

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(tmp), **kwargs)

        def log_message(self, fmt: str, *args: object) -> None:
            return

    socketserver.TCPServer.allow_reuse_address = True
    httpd = socketserver.TCPServer(("127.0.0.1", 0), Handler)
    httpd.timeout = 2
    port = httpd.server_address[1]
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{port}/scfm-reg-audit"
    try:
        paths = [
            "/",
            "/figures/",
            "/reproducibility/",
        ]
        for path in paths:
            url = base + path
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=5) as resp:
                code = resp.status
                if code != 200:
                    fail(f"smoke {url} -> {code}")
            ok(f"smoke 200 {path}")
        missing = base + "/no-such-path/"
        try:
            urllib.request.urlopen(missing, timeout=5)
        except urllib.error.HTTPError as err:
            if err.code not in {404, 403}:
                fail(f"smoke missing path -> {err.code}")
            ok(f"smoke missing path {err.code}")
    finally:
        httpd.shutdown()
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
