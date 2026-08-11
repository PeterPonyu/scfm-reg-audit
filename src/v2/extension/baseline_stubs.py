#!/usr/bin/env python3
"""Simple baseline method stubs for extension Tier A–C comparisons.

Emits dry-run plans only. No FM BH membership; no PeerJ Support writes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from registry import HEAVY_ARTIFACT_ROOT, load_extension_registry  # noqa: E402

ROOT = _EXT.parents[2]
PLAN_ROOT = ROOT / "docs" / "reports" / "extension-plans" / "baselines"

STUB_METHODS = (
    "degree_matched_random",
    "motif_only_rp",
    "encode_chip_binding",
    "collectri_prior",
)


def baseline_plan(method_id: str) -> dict[str, Any]:
    reg = load_extension_registry()
    meta = reg.get_method(method_id)
    if meta.get("bh_membership") == "fm_families":
        raise PermissionError(f"{method_id} is an FM-family method, not a baseline stub")
    return {
        "method_id": method_id,
        "tier": meta.get("tier"),
        "estimand": meta.get("estimand"),
        "bh_membership": meta.get("bh_membership"),
        "status": meta.get("status", "stub"),
        "panel": "frozen_446x1200",
        "outputs": {
            "score_matrix": f"{HEAVY_ARTIFACT_ROOT}baselines/{method_id}/scores.npz",
            "summary": f"{HEAVY_ARTIFACT_ROOT}baselines/{method_id}/summary.json",
        },
        "peerj_bh_families_unchanged": 8,
        "hooks": {
            "degree_matched_random": "sample edges with degree-matched null on panel mask",
            "motif_only_rp": "reuse motif hits; ablate ATAC accessibility weights",
            "encode_chip_binding": "project local ENCODE peaks onto panel; extend 0x40 SI",
            "collectri_prior": "load cached CollecTRI/OmniPath if present; else skip",
        }.get(method_id, "implement score matrix on panel"),
        "notes": meta.get("notes"),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=STUB_METHODS, default="motif_only_rp")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)

    methods = list(STUB_METHODS) if args.all else [args.method]
    plans = [baseline_plan(mid) for mid in methods]
    print(json.dumps(plans if args.all else plans[0], indent=2))

    if args.write:
        for plan in plans:
            out = PLAN_ROOT / plan["method_id"]
            out.mkdir(parents=True, exist_ok=True)
            path = out / "baseline_plan.dry_run.json"
            path.write_text(json.dumps(plan, indent=2) + "\n")
            print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
