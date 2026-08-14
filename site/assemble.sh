#!/usr/bin/env bash
# Build the Hugo site and copy the two SoT PDFs into the publish dir (artifact only).
# Fail closed if a PDF is missing. Do not copy FigureN.pdf. Do not link them from pages.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE="${ROOT}/site"
OUT="${ASSEMBLE_OUT:-${SITE}/public}"
HUGO_BIN="${HUGO:-hugo}"

die() {
  printf 'assemble: %s\n' "$*" >&2
  exit 1
}

command -v "${HUGO_BIN}" >/dev/null 2>&1 || die "hugo not on PATH (set HUGO=)"
ver="$("${HUGO_BIN}" version)"
printf '%s\n' "${ver}" | grep -q 'extended' || die "need Hugo extended, got: ${ver}"

PEERJ_SRC="${ROOT}/paper/submission_peerj/flat_upload/manuscript.pdf"
FRONTIERS_SRC="${ROOT}/paper/submission_frontiers_genetics/manuscript.pdf"
[[ -f "${PEERJ_SRC}" ]] || die "missing ${PEERJ_SRC}"
[[ -f "${FRONTIERS_SRC}" ]] || die "missing ${FRONTIERS_SRC}"
[[ -f "${SITE}/hugo.toml" ]] || die "missing ${SITE}/hugo.toml"

rm -rf "${OUT}"
mkdir -p "${OUT}"

(
  cd "${SITE}"
  "${HUGO_BIN}" --minify --destination "${OUT}"
)

python3 - "${OUT}" "${PEERJ_SRC}" "${FRONTIERS_SRC}" <<'PY'
import shutil
import sys
from pathlib import Path

out = Path(sys.argv[1])
peerj_src = Path(sys.argv[2])
frontiers_src = Path(sys.argv[3])
dest_dir = out / "products"
dest_dir.mkdir(parents=True, exist_ok=True)
peerj = dest_dir / "peerj-manuscript.pdf"
frontiers = dest_dir / "frontiers-manuscript.pdf"
shutil.copy2(peerj_src, peerj)
shutil.copy2(frontiers_src, frontiers)
if peerj.stat().st_size != peerj_src.stat().st_size:
    sys.exit("assemble: PeerJ PDF size mismatch after copy")
if frontiers.stat().st_size != frontiers_src.stat().st_size:
    sys.exit("assemble: Frontiers PDF size mismatch after copy")

pdfs = sorted(p for p in out.rglob("*") if p.is_file() and p.suffix.lower() == ".pdf")
if len(pdfs) != 2:
    sys.exit(f"assemble: expected exactly 2 PDFs, found {len(pdfs)}")
names = {p.name for p in pdfs}
if names != {"peerj-manuscript.pdf", "frontiers-manuscript.pdf"}:
    sys.exit(f"assemble: unexpected PDF names {sorted(names)}")

deny_sub = (
    "docs/reports",
    "results/v2",
    "visual_qa",
    "flat_upload",
    ".omx/",
    ".omc/",
)
deny_suf = (".h5ad", ".h5", ".npz")
deny_names = tuple(f"Figure{i}.pdf" for i in range(1, 14)) + (
    "FigureA1.pdf",
    "FigureA2.pdf",
    "FigureA3.pdf",
)
for path in out.rglob("*"):
    if not path.is_file():
        continue
    rel = path.relative_to(out).as_posix()
    if any(n in rel for n in deny_sub) or rel.endswith(deny_suf) or path.name in deny_names:
        sys.exit(f"assemble: denylist hit {rel!r}")

print(f"assemble output: {out}")
files = [p for p in out.rglob("*") if p.is_file()]
print(f"assemble files: {len(files)}")
for path in sorted(files):
    print(f"  {path.relative_to(out).as_posix()}")
PY
