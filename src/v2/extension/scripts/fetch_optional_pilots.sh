#!/usr/bin/env bash
# Optional small cancer/dev pilot fetches — ESTIMATE / DRY-RUN by default.
# This session constraint: do NOT download multi-GB atlases.
#
# Usage:
#   ./src/v2/extension/scripts/fetch_optional_pilots.sh
#   CONFIRM_FETCH=1 DESCARTES_TISSUE_RDS_URL=... ./src/v2/extension/scripts/...
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
PILOT_ROOT="${EXTENSION_PILOT_ROOT:-${DESKTOP_DATA:-${HOME}/Desktop/data}/datasets/extension_pilots}"
MODE="${1:-estimate}"

echo "[extension-fetch] mode=$MODE"
echo "[extension-fetch] repo=$ROOT"
echo "[extension-fetch] pilot_root set via EXTENSION_PILOT_ROOT or DESKTOP_DATA"
echo
echo "=== Order-of-magnitude size estimates (no download) ==="
cat <<'EOF'
| Scenario                         | Network     | Disk after   | Notes |
|----------------------------------|------------:|-------------:|-------|
| Tiny construct SI (DESCARTES spleen RDS) | ~0.1 GB | ~0.1–0.2 GB | Preferred construct SI; URL must be set |
| BMMC multiome (already local)    | 0 GB        | ~5.7 GB      | Compute-only; panel policy gated |
| Modest HTAN open sample pilot    | ~0.5–5 GB   | ~1–6 GB      | One sample / few matrices; Synapse open only |
| Whole HTAN / DESCARTES RAW lake  | 10s–100s GB | same         | FORBIDDEN without explicit go-ahead |
| Cancer/Dev RNA lakes 27/28       | 0 (local)   | ~90 GB       | out_of_scope for G_ATAC / Support |
EOF
echo
echo "See docs/reports/optional_cancer_dev_download_costs.md for 5090 VRAM/time."

if [[ "${CONFIRM_FETCH:-}" != "1" ]]; then
  echo
  echo "[dry-run] CONFIRM_FETCH!=1 → not downloading anything."
  echo "[dry-run] Set DESCARTES_TISSUE_RDS_URL + CONFIRM_FETCH=1 to pull one ~100MB RDS."
  exit 0
fi

if [[ -z "${DESCARTES_TISSUE_RDS_URL:-}" ]]; then
  echo "ABORT: CONFIRM_FETCH=1 requires DESCARTES_TISSUE_RDS_URL" >&2
  exit 2
fi

mkdir -p "$PILOT_ROOT/descartes"
dest="$PILOT_ROOT/descartes/tissue_pilot.seurat.RDS.gz"
echo "[fetch] $DESCARTES_TISSUE_RDS_URL → $dest"
curl -fL --continue-at - -o "$dest" "$DESCARTES_TISSUE_RDS_URL"
sz=$(stat -c%s "$dest")
# Hard stop: >2GB unexpected for "tiny" pilot
if (( sz > 2000000000 )); then
  echo "ABORT: downloaded file >2GB ($sz bytes); removing" >&2
  rm -f "$dest"
  exit 3
fi
echo "[ok] fetched $sz bytes"
