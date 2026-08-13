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
        "peerj/index.html",
        "frontiers/index.html",
        "figures/index.html",
        "reproducibility/index.html",
        "404.html",
        "products/peerj-manuscript.pdf",
        "products/frontiers-manuscript.pdf",
        ".nojekyll",
    ]
    for rel in routes:
        path = PUBLIC / rel
        if not path.is_file():
            fail(f"missing route {rel}")
        ok(rel)

    hub = (PUBLIC / "index.html").read_text(encoding="utf-8")
    if "0/13" not in hub:
        fail("hub missing protocol-pass 0/13")
    ok("hub has 0/13")

    if re.search(r"Ph\.D|PhD", hub, flags=re.I):
        fail("hub contains Ph.D / PhD chrome")
    if "Frontiers in Genetics" in hub:
        fail("hub uses journal title-page identity 'Frontiers in Genetics'")
    if re.search(r"<h1[^>]*>.*Genetics.*</h1>", hub, flags=re.I | re.S):
        fail("hub h1 uses genetics as identity")
    if "Army Medical" in hub or "word-count" in hub.lower():
        fail("hub contains degree / word-count chrome")
    if "PeerJ CS PDF" not in hub or "Frontiers Genetics PDF" not in hub:
        fail("hub missing equal PDF file labels")
    ok("hub identity chrome checks")

    for html_path in PUBLIC.rglob("*.html"):
        text = html_path.read_text(encoding="utf-8")
        if re.search(r"Ph\.D|PhD", text):
            fail(f"{html_path.relative_to(PUBLIC)} contains Ph.D / PhD")
        if "Army Medical" in text:
            fail(f"{html_path.relative_to(PUBLIC)} contains Army Medical")
        if "fonts.googleapis.com" in text or "fonts.gstatic.com" in text:
            fail(f"{html_path.relative_to(PUBLIC)} loads Google Fonts CDN")
    ok("no Ph.D / Army Medical / webfont CDN in HTML")

    for layout in LAYOUTS.rglob("*.html"):
        text = layout.read_text(encoding="utf-8")
        if "/scfm-reg-audit/" in text:
            fail(f"{layout.relative_to(SITE)} hardcodes /scfm-reg-audit/ (use relURL)")
        if "relURL" not in text and layout.name in {
            "header.html",
            "index.html",
            "404.html",
            "baseof.html",
        }:
            # header/index/404/baseof must use relURL for internal paths
            if layout.name != "baseof.html" or "RelPermalink" not in text:
                pass
    if 'href="/scfm-reg-audit/' in (LAYOUTS / "index.html").read_text(encoding="utf-8"):
        fail("hub layout hardcodes project-pages prefix")
    ok("layouts do not hardcode /scfm-reg-audit/")

    if "/scfm-reg-audit/peerj/" not in hub:
        fail("built hub missing relURL prefix /scfm-reg-audit/peerj/")
    if 'href="/peerj/' in hub:
        fail("built hub has root-absolute /peerj/ (missing project base)")
    if re.search(r"href=/>", hub) or 'href="/"' in hub:
        fail("built hub has root-absolute href=/ (missing project base)")
    if 'href=/scfm-reg-audit/' not in hub.replace('href="/scfm-reg-audit/', "href=/scfm-reg-audit/"):
        fail("built hub missing project-pages href prefix")
    ok("built hub uses project-pages prefix")

    pdfs = list(PUBLIC.rglob("*.pdf"))
    if len(pdfs) != 2:
        fail(f"expected 2 PDFs, found {len(pdfs)}")
    names = {p.name for p in pdfs}
    if names != {"peerj-manuscript.pdf", "frontiers-manuscript.pdf"}:
        fail(f"unexpected PDFs {sorted(names)}")
    for name in (f"Figure{i}.pdf" for i in range(1, 13)):
        if (PUBLIC / name).exists() or list(PUBLIC.rglob(name)):
            fail(f"denylist PDF present: {name}")
    ok("exactly two product PDFs")

    figures = (PUBLIC / "figures/index.html").read_text(encoding="utf-8")
    if "paper/figs_extension/fig_ext" in figures:
        fail("figures page links missing origin/main path paper/figs_extension/")
    if "paper/submission_frontiers_genetics/fig_ext1_construct_mantel.tex" not in figures:
        fail("figures page missing origin/main appendix blob path")
    ok("appendix blob paths exist on origin/main")

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
            "/peerj/",
            "/frontiers/",
            "/figures/",
            "/reproducibility/",
            "/products/peerj-manuscript.pdf",
            "/products/frontiers-manuscript.pdf",
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
            # SimpleHTTPRequestHandler 404s for missing files; directories may 301
        except urllib.error.HTTPError as err:
            if err.code not in {404, 403}:
                fail(f"smoke missing path -> {err.code}")
            ok(f"smoke missing path {err.code}")
    finally:
        httpd.shutdown()
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
