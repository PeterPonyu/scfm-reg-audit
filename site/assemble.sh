#!/usr/bin/env bash
# Assemble the GitHub Pages artifact from site/ source + two product PDFs.
# Output: site/_site/ (gitignored). Fail closed if a SoT PDF is missing.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE="${ROOT}/site"
OUT="${ASSEMBLE_OUT:-${SITE}/_site}"
ALLOWLIST="${SITE}/allowlist.json"

die() {
  printf 'assemble: %s\n' "$*" >&2
  exit 1
}

[[ -f "${ALLOWLIST}" ]] || die "missing ${ALLOWLIST}"

PEERJ_SRC="${ROOT}/paper/manuscript.pdf"
FRONTIERS_SRC="${ROOT}/paper/submission_frontiers_genetics/manuscript.pdf"
[[ -f "${PEERJ_SRC}" ]] || die "missing ${PEERJ_SRC}"
[[ -f "${FRONTIERS_SRC}" ]] || die "missing ${FRONTIERS_SRC}"

rm -rf "${OUT}"
mkdir -p "${OUT}"

python3 - "${SITE}" "${OUT}" "${ROOT}" "${ALLOWLIST}" "${PEERJ_SRC}" "${FRONTIERS_SRC}" <<'PY'
import json
import shutil
import sys
from pathlib import Path

site = Path(sys.argv[1])
out = Path(sys.argv[2])
root = Path(sys.argv[3])
allow = json.loads(Path(sys.argv[4]).read_text())
peerj_src = Path(sys.argv[5])
frontiers_src = Path(sys.argv[6])

for rel in allow["copy"]:
    src = site / rel
    if not src.is_file():
        sys.exit(f"assemble: missing source {rel}")
    dest = out / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)

pdfs = allow["pdfs"]
if len(pdfs) != 2:
    sys.exit("assemble: allowlist must name exactly two PDFs")
expected = [
    ("paper/manuscript.pdf", "products/peerj-manuscript.pdf"),
    ("paper/submission_frontiers_genetics/manuscript.pdf", "products/frontiers-manuscript.pdf"),
]
got = [(item["src"], item["dest"]) for item in pdfs]
if got != expected:
    sys.exit(f"assemble: PDF copy-set mismatch: {got}")

dest_peerj = out / "products/peerj-manuscript.pdf"
dest_frontiers = out / "products/frontiers-manuscript.pdf"
dest_peerj.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(peerj_src, dest_peerj)
shutil.copy2(frontiers_src, dest_frontiers)
if dest_peerj.stat().st_size != (root / "paper/manuscript.pdf").stat().st_size:
    sys.exit("assemble: PeerJ PDF size mismatch after copy")
if dest_frontiers.stat().st_size != (
    root / "paper/submission_frontiers_genetics/manuscript.pdf"
).stat().st_size:
    sys.exit("assemble: Frontiers PDF size mismatch after copy")

pdf_files = sorted(p for p in out.rglob("*") if p.is_file() and p.suffix.lower() == ".pdf")
if len(pdf_files) != 2:
    sys.exit(f"assemble: expected exactly 2 PDFs in artifact, found {len(pdf_files)}")
names = {p.name for p in pdf_files}
if names != {"peerj-manuscript.pdf", "frontiers-manuscript.pdf"}:
    sys.exit(f"assemble: unexpected PDF names {sorted(names)}")

deny_sub = allow["deny_substrings"]
deny_suf = tuple(allow["deny_suffixes"])
for path in out.rglob("*"):
    if not path.is_file():
        continue
    rel = path.relative_to(out).as_posix()
    for needle in deny_sub:
        if needle in rel:
            sys.exit(f"assemble: denylist hit {rel!r} ({needle})")
    if rel.endswith(deny_suf):
        sys.exit(f"assemble: denylist suffix {rel!r}")

print(f"assemble output: {out}")
files = [p for p in out.rglob("*") if p.is_file()]
print(f"assemble files: {len(files)}")
for path in sorted(files):
    print(f"  {path.relative_to(out).as_posix()}")
PY
