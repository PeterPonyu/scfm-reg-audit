#!/usr/bin/env bash
# Assemble then run publish-dir checks (science routes, 0/13, banned-word invert, smoke).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "${ROOT}/site/assemble.sh"
python3 "${ROOT}/site/test.py" "${ROOT}/site/public"
