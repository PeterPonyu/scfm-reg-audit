#!/usr/bin/env bash
# Thin wrapper around fetch_approved_plan.py (approval-gated).
# Requires: SCREG_DOWNLOAD_APPROVED=1 and SCREG_DOWNLOAD_PLAN_ID=<plan>
# Default: dry-run. Pass --execute only after verbal human approval.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
PLAN_ID="${1:?usage: fetch_plan.sh <D4|D5|...> [--execute] [--asset-id ID]}"
shift || true
exec python "${ROOT}/src/v2/extension/scripts/fetch_approved_plan.py" \
  --plan-id "${PLAN_ID}" "$@"
