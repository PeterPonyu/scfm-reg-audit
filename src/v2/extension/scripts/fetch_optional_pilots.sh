#!/usr/bin/env bash
# DEPRECATED as a download constructor (Option B-prime / PeerJ freeze wave).
# Downloads are approval-document-only. See:
#   docs/reports/download_approval_optional_pilots.md
#   docs/reports/optional_cancer_dev_download_costs.md
#
# This script intentionally refuses to fetch. No network downloader path ships here.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
APPROVAL="${ROOT}/docs/reports/download_approval_optional_pilots.md"
COSTS="${ROOT}/docs/reports/optional_cancer_dev_download_costs.md"

echo "[extension-fetch] DOWNLOAD CONSTRUCTOR DISABLED"
echo "[extension-fetch] Approval document: ${APPROVAL}"
echo "[extension-fetch] Cost / policy model: ${COSTS}"
echo
echo "No network fetch will be performed from this script."
echo "After human approval, fetch manually outside the PeerJ freeze path;"
echo "construct code lives in src/v2/extension/ (local G_ATAC assets only)."
exit 0
