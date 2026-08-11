#!/usr/bin/env python3
"""Fail-closed download seam for extension pilots (no network I/O).

Requires ``SCREG_DOWNLOAD_APPROVED=1`` and a matching plan id. Even when
approved, this module never fetches; it only prints a manual recipe or writes
a dry-run plan JSON under ``results/v2/extension/download-plans/``.

The env gate is an **operator checklist ceremony**, not a security or auth
boundary: any local process can set the same variables. Do not add a real
network fetch behind these env vars without a stronger binder (e.g. signed
approval artifact, human-in-the-loop out-of-band confirm, or capability that
cannot be forged by env alone).

Exit code **0** means “recipe emitted / dry-run only” — never “download
succeeded”. Exit code **2** means the gate refused or the plan id was invalid.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from paths import (  # noqa: E402
    EXTENSION_ROOT,
    ROOT,
    assert_confined_write_path,
    assert_safe_tag,
)

APPROVAL_ENV = "SCREG_DOWNLOAD_APPROVED"
PLAN_MATCH_ENV = "SCREG_DOWNLOAD_PLAN_ID"
APPROVAL_DOC = ROOT / "docs" / "reports" / "download_approval_optional_pilots.md"
COSTS_DOC = ROOT / "docs" / "reports" / "optional_cancer_dev_download_costs.md"

# Keep in sync with the approval matrix in download_approval_optional_pilots.md.
PLAN_REGISTRY: dict[str, dict[str, Any]] = {
    "D0": {
        "title": "Nothing (local-only)",
        "decision": "default",
        "fetchable": False,
        "rejected": False,
        "manual_recipe": [
            "No network fetch — use existing local fibro/brain/PBMC G_ATAC under results/v2/.",
            "Construct / baselines / claim-pack only via src/v2/extension/cli.py.",
        ],
    },
    "D1": {
        "title": "DESCARTES spleen RDS (tiny construct SI)",
        "decision": "pending",
        "fetchable": True,
        "rejected": False,
        "manual_recipe": [
            "Record exact GEO/URL in the approval checklist.",
            "Fetch manually (outside this CLI) into ${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/.",
            "Expected compressed size gate: < 2 GB (target ~0.1 GB).",
            "Convert RDS → peak h5ad (var_names chrom:start-end) as descartes_spleen_peaks.h5ad in that dir.",
            "Validate: python src/v2/extension/cli.py descartes-bridge",
            "Build overlay: TAG=DESCARTES_spleen ATAC_FILE=<h5ad> SCREG_EXTENSION_OUT=results/v2/extension/construct/DESCARTES_spleen SCREG_PEERJ_SUPPORT_LOCK=1 python src/v2/build_atac_graph_v2.py",
            "Then: python src/v2/extension/cli.py construct --tissue descartes_spleen --execute --write",
            "Never write into PeerJ Support JSON / canonical results/v2/ G_ATAC locks.",
        ],
    },
    "D2": {
        "title": "BMMC multiome (already local)",
        "decision": "no_new_download",
        "fetchable": False,
        "rejected": False,
        "manual_recipe": [
            "No new download — BMMC h5ad should already be under ${DESKTOP_DATA}/external/scfm-reg-audit/gse194122/.",
            "Set SCREG_BMMC_H5AD if needed; P3 = construct only (no FM Support).",
            "Prepare peaks: python src/v2/extension/cli.py prepare-bmmc --execute --write",
            "Build overlay G_ATAC via the build_command emitted by prepare-bmmc (SCREG_EXTENSION_OUT + SCREG_PEERJ_SUPPORT_LOCK=1).",
            "Then: python src/v2/extension/cli.py construct --tissue bmmc --execute --write",
        ],
    },
    "D3": {
        "title": "HTAN open single-sample pilot",
        "decision": "pending",
        "fetchable": True,
        "rejected": False,
        "manual_recipe": [
            "Confirm Synapse/open access + exact sample ID in the approval checklist.",
            "Fetch manually into ${DESKTOP_DATA}/datasets/extension_pilots/htan_pilot/.",
            "Size gate: confirm < 10 GB compressed before fetch.",
            "Extension-only artifacts under results/v2/extension/; never inflate PeerJ 13-row SAP.",
        ],
    },
    "D4": {
        "title": "Whole HTAN / DESCARTES RAW / File_S6 lakes",
        "decision": "rejected",
        "fetchable": False,
        "rejected": True,
        "manual_recipe": ["Rejected — must refuse lake-scale fetch for Support / G_ATAC."],
    },
    "D5": {
        "title": "Cancer 28 / Dev 27 RNA lakes",
        "decision": "rejected",
        "fetchable": False,
        "rejected": True,
        "manual_recipe": [
            "Rejected — RNA∩ATAC empty; forbidden for Support / G_ATAC."
        ],
    },
}


def known_plan_ids() -> frozenset[str]:
    return frozenset(PLAN_REGISTRY)


def _approval_granted() -> bool:
    return os.environ.get(APPROVAL_ENV, "").strip() == "1"


def _matched_plan_id(plan_id: str) -> bool:
    env_plan = os.environ.get(PLAN_MATCH_ENV, "").strip()
    return bool(env_plan) and env_plan == plan_id


def build_dry_run_plan(plan_id: str) -> dict[str, Any]:
    meta = PLAN_REGISTRY[plan_id]
    return {
        "schema_version": 1,
        "plan_id": plan_id,
        "title": meta["title"],
        "decision": meta["decision"],
        "fetchable": meta["fetchable"],
        "rejected": meta["rejected"],
        "network_fetch_performed": False,
        "approval_doc": str(APPROVAL_DOC.relative_to(ROOT)),
        "costs_doc": str(COSTS_DOC.relative_to(ROOT)),
        "manual_recipe": list(meta["manual_recipe"]),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "note": (
            "Fail-closed stub: no urllib/wget/curl. After human approval + env gates, "
            "execute the manual_recipe outside this CLI."
        ),
    }


def run_download_gate(plan_id: str, *, write_plan: bool = True) -> tuple[int, dict[str, Any]]:
    """Gate download by env + plan id. Never performs network I/O.

    Returns (exit_code, payload). Exit 0 = dry-run recipe only; exit 2 = refused.
    """
    try:
        plan_id = assert_safe_tag(plan_id, label="plan-id")
    except ValueError as exc:
        payload = {
            "status": "invalid_plan_id",
            "plan_id": plan_id,
            "message": str(exc),
            "network_fetch_performed": False,
        }
        return 2, payload

    if plan_id not in PLAN_REGISTRY:
        payload = {
            "status": "unknown_plan",
            "plan_id": plan_id,
            "known_plan_ids": sorted(PLAN_REGISTRY),
            "message": (
                f"Unknown plan id {plan_id!r}. See {APPROVAL_DOC.relative_to(ROOT)} "
                f"(known: {', '.join(sorted(PLAN_REGISTRY))})."
            ),
            "network_fetch_performed": False,
        }
        return 2, payload

    meta = PLAN_REGISTRY[plan_id]
    if meta["rejected"]:
        payload = {
            "status": "rejected_plan",
            "plan_id": plan_id,
            "message": (
                f"Plan {plan_id} is permanently rejected in the approval matrix. "
                "Refuse lake-scale / RNA-lake fetch."
            ),
            "network_fetch_performed": False,
            "manual_recipe": list(meta["manual_recipe"]),
        }
        return 2, payload

    if not _approval_granted() or not _matched_plan_id(plan_id):
        payload = {
            "status": "approval_required",
            "plan_id": plan_id,
            "message": (
                "Download gate closed. Set SCREG_DOWNLOAD_APPROVED=1 and "
                f"SCREG_DOWNLOAD_PLAN_ID={plan_id} after filling "
                f"{APPROVAL_DOC.relative_to(ROOT)}. This CLI never fetches."
            ),
            "required_env": {
                APPROVAL_ENV: "1",
                PLAN_MATCH_ENV: plan_id,
            },
            "network_fetch_performed": False,
        }
        return 2, payload

    plan = build_dry_run_plan(plan_id)
    plan["status"] = "approved_dry_run_no_fetch"
    if write_plan:
        out_dir = assert_confined_write_path(
            EXTENSION_ROOT / "download-plans",
            label="download-plans",
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        plan_path = out_dir / f"{plan_id}.dry_run.json"
        plan_path.write_text(json.dumps(plan, indent=2) + "\n")
        plan["dry_run_plan_path"] = str(plan_path.relative_to(ROOT))

    plan["next_manual_commands"] = [
        f"# Approved dry-run only — execute manually if still desired:",
        *[f"#   {step}" for step in plan["manual_recipe"]],
        f"# Approval docs: {APPROVAL_DOC.relative_to(ROOT)}",
    ]
    return 0, plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fail-closed download gate (no network I/O). "
            "Exit 0 = recipe emitted / dry-run only — never means download succeeded. "
            "Env gates (SCREG_DOWNLOAD_*) are an operator checklist ceremony, "
            "not a security/auth boundary."
        )
    )
    parser.add_argument(
        "--plan-id",
        required=True,
        help="Approval matrix id (D0–D5) from download_approval_optional_pilots.md",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Do not write dry-run JSON under results/v2/extension/download-plans/",
    )
    args = parser.parse_args(argv)
    code, payload = run_download_gate(args.plan_id, write_plan=not args.no_write)
    print(json.dumps(payload, indent=2))
    if code != 0:
        print(payload.get("message", "download gate closed"), file=sys.stderr)
    else:
        print(
            "Exit 0 = recipe emitted / dry-run only — never means download succeeded. "
            "No network fetch performed.",
            file=sys.stderr,
        )
        for line in payload.get("next_manual_commands") or []:
            print(line, file=sys.stderr)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
