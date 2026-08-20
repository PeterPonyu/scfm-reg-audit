#!/usr/bin/env bash
# Build the Hugo site. Do not copy FigureN.pdf or full manuscripts into the publish dir.
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

[[ -f "${SITE}/hugo.toml" ]] || die "missing ${SITE}/hugo.toml"

rm -rf "${OUT}"
mkdir -p "${OUT}"

(
  cd "${SITE}"
  "${HUGO_BIN}" --minify --destination "${OUT}"
)

python3 - "${OUT}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
pdfs = sorted(p for p in out.rglob("*") if p.is_file() and p.suffix.lower() == ".pdf")
if pdfs:
    names = sorted(p.relative_to(out).as_posix() for p in pdfs)
    sys.exit(f"assemble: expected no PDFs in the publish dir, found {names}")

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
