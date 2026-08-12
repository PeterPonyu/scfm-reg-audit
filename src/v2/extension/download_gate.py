#!/usr/bin/env python3
"""Fail-closed download seam for extension pilots (no network I/O).

Requires ``SCREG_DOWNLOAD_APPROVED=1`` and a matching plan id. Even when
approved, this module never fetches; it only prints a manual recipe or writes
a dry-run plan JSON under ``results/v2/extension/download-plans/``.

Real network fetch (when later approved) lives in
``src/v2/extension/scripts/fetch_approved_plan.py`` — never in this CLI path.

The env gate is an **operator checklist ceremony**, not a security or auth
boundary: any local process can set the same variables. Do not add a real
network fetch behind these env vars without a stronger binder (e.g. signed
approval artifact, human-in-the-loop out-of-band confirm, or capability that
cannot be forged by env alone).

Exit code **0** means “recipe emitted / dry-run only” — never “download
succeeded”. Exit code **2** means the gate refused or the plan id was invalid.

D4/D5 are **pending_large** (approval-gated infrastructure), not permanently
impossible. Prior “cannot download” language was **policy** (PeerJ freeze /
disk / G_ATAC identity), not a technical block.
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
FETCH_LOG_HINT = (
    "${DESKTOP_DATA}/datasets/extension_pilots/manifests/FETCH_LOG.md"
)

# Keep in sync with the approval matrix in download_approval_optional_pilots.md.
PLAN_REGISTRY: dict[str, dict[str, Any]] = {
    "D0": {
        "title": "Nothing (local-only)",
        "decision": "default",
        "fetchable": False,
        "rejected": False,
        "assets": [],
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
        "assets": [
            {
                "id": "GSM4508940_spleen_rds_gz",
                "url": (
                    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4508nnn/"
                    "GSM4508940/suppl/GSM4508940_spleen_filtered.seurat.RDS.gz"
                ),
                "dest": (
                    "${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/"
                    "GSM4508940_spleen_filtered.seurat.RDS.gz"
                ),
                "bytes_estimate": 99_397_382,
                "size_label": "~0.1 GB",
                "checksum_sha256": (
                    "c024a30868e511e1dc7b2e72efe1700bf7033a2af85ee2c341dd4ee1b950e397"
                ),
                "checksum_note": "Recorded in FETCH_LOG after 2026-08-12 D1 fetch.",
            }
        ],
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
        "assets": [],
        "manual_recipe": [
            "No new download — BMMC h5ad should already be under ${DESKTOP_DATA}/external/scfm-reg-audit/gse194122/.",
            "Set SCREG_BMMC_H5AD if needed; P3 = construct only (no FM Support).",
            "Prepare peaks: python src/v2/extension/cli.py prepare-bmmc --execute --write",
            "Build overlay G_ATAC via the build_command emitted by prepare-bmmc (SCREG_EXTENSION_OUT + SCREG_PEERJ_SUPPORT_LOCK=1).",
            "Then: python src/v2/extension/cli.py construct --tissue bmmc --execute --write",
        ],
    },
    "D3": {
        "title": "HTAN open single-sample pilot (local tar)",
        "decision": "no_new_download",
        "fetchable": False,
        "rejected": False,
        "assets": [],
        "manual_recipe": [
            "No new network download — GSM7710026 C3N-01334 snATAC tar is already local.",
            "Canonical root: ${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/",
            "Inventory / prepare (not download): python src/v2/extension/cli.py prepare-htan --write",
            "Tar is Cell Ranger fragments-only → structured status=blocked (E4.3b) unless a peak h5ad is supplied.",
            "If peak h5ad appears with chrom:start-end var_names, re-run prepare-htan then emitted build_command + construct --tissue htan_gbm_pilot.",
            "Extension-only artifacts under results/v2/extension/; never inflate PeerJ 13-row SAP.",
            "Do not implement this path via download --plan-id D3 execute/fetch.",
        ],
    },
    "D4": {
        "title": "Whole HTAN / DESCARTES RAW / File_S6 lakes",
        "decision": "pending_large",
        "fetchable": True,
        "rejected": False,
        "policy_risk": [
            "Unbounded lake / exceeds pilot hard gates (was policy-rejected, not technically impossible).",
            "PeerJ freeze narrative + disk risk (~29 GB RAW + ~4.3 GB File_S6; HTAN cohort TB-scale).",
            "Must not inflate PeerJ Support / locked G_ATAC; extension overlay only after subsetting.",
        ],
        "assets": [
            {
                "id": "GSE149683_RAW",
                "url": (
                    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/"
                    "GSE149683/suppl/GSE149683_RAW.tar"
                ),
                "dest": (
                    "${DESKTOP_DATA}/datasets/extension_pilots/descartes_lake/"
                    "GSE149683_RAW.tar"
                ),
                "bytes_estimate": 31_177_328_640,
                "size_label": "~29.0 GB",
                "checksum_note": (
                    "No published sha256 in FETCH_LOG; verify size == 31177328640 "
                    "(GEO filelist / descartes/filelist.txt) after fetch."
                ),
            },
            {
                "id": "GSE149683_File_S6",
                "url": (
                    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/"
                    "GSE149683/suppl/"
                    "GSE149683_File_S6.Cicero_gene_activity_scores_by_cell_type.csv.gz"
                ),
                "dest": (
                    "${DESKTOP_DATA}/datasets/extension_pilots/descartes_lake/"
                    "GSE149683_File_S6.Cicero_gene_activity_scores_by_cell_type.csv.gz"
                ),
                "bytes_estimate": 4_618_000_000,
                "size_label": "~4.3 GB",
                "checksum_note": (
                    "GEO lists ~4.3 Gb; no sha256 in FETCH_LOG — record sha256 after fetch."
                ),
            },
            {
                "id": "HTAN_cohort_wide",
                "url": None,
                "dest": (
                    "${DESKTOP_DATA}/datasets/extension_pilots/htan_lake/"
                ),
                "bytes_estimate": None,
                "size_label": "TB-scale (Synapse cohort)",
                "checksum_note": (
                    "No single HTTP URL. Shortlist open synIDs from "
                    "${DESKTOP_DATA}/datasets/extension_pilots/htan/"
                    "Sample_ID_Lookup_table_in_repositories.xlsx then "
                    "`synapse get <synID>` — never auto-dump cohort-wide."
                ),
                "fetch_via": "synapse_manual",
            },
        ],
        "manual_recipe": [
            "POLICY: previously blocked as unbounded lake (not technical impossibility).",
            "Dual gate: (1) verbal approve D4 + SCREG_DOWNLOAD_* env; (2) disk ≥ ~40 GB free for RAW+S6.",
            "Dry-run recipe only via CLI; real fetch via scripts/fetch_approved_plan.py --execute.",
            "Destinations under ${DESKTOP_DATA}/datasets/extension_pilots/descartes_lake/ (not PeerJ package).",
            "Optional File_S5 (~392 MB) is NOT in default D4 asset list — add only if explicitly requested.",
            "HTAN whole-lake is Synapse/TB — not emitted as an auto-URL; keep to shortlisted open samples.",
            "After fetch: construct-only / lake_blocked registry roles; never PeerJ Support / G_ATAC identity.",
            f"Log bytes+sha256 into {FETCH_LOG_HINT}.",
        ],
    },
    "D5": {
        "title": "Cancer 28 / Dev 27 RNA lakes",
        "decision": "pending_large",
        "fetchable": True,
        "rejected": False,
        "policy_risk": [
            "RNA∩ATAC empty — forbidden for Support / G_ATAC identity even after download.",
            "Largely already local (~90 GB class); re-fetch usually pointless.",
            "PeerJ freeze: do not wire lakes into capsule validation or Support JSON.",
        ],
        "assets": [
            {
                "id": "cancer28_local",
                "url": None,
                "dest": "${DESKTOP_DATA}/datasets/CancerDatasets/",
                "bytes_estimate": 60_000_000_000,
                "size_label": "~60 GB already local (28 h5ads + CancerDatasets2)",
                "local_dirs": [
                    "${DESKTOP_DATA}/datasets/CancerDatasets",
                    "${DESKTOP_DATA}/datasets/CancerDatasets2",
                ],
                "checksum_note": "Inventory-only; no GEO re-fetch by default.",
                "fetch_via": "local_inventory",
            },
            {
                "id": "dev27_local",
                "url": None,
                "dest": "${DESKTOP_DATA}/datasets/DevelopmentDatasets/",
                "bytes_estimate": 31_000_000_000,
                "size_label": "~31 GB already local (27 h5ads + DevelopmentDatasets2)",
                "local_dirs": [
                    "${DESKTOP_DATA}/datasets/DevelopmentDatasets",
                    "${DESKTOP_DATA}/datasets/DevelopmentDatasets2",
                ],
                "checksum_note": "Inventory-only; no GEO re-fetch by default.",
                "fetch_via": "local_inventory",
            },
        ],
        "manual_recipe": [
            "POLICY: blocked from Support/G_ATAC (RNA-only estimand mismatch), not because files cannot exist locally.",
            "Default action after approve: inventory local lakes (no network).",
            "Paths: ${DESKTOP_DATA}/datasets/CancerDatasets{,2}/ and DevelopmentDatasets{,2}/.",
            "EVEN AFTER DOWNLOAD/INVENTORY: must not enter Support / G_ATAC construction.",
            "Coexpression-side / Limitations use only; registry roles remain out_of_scope / rna_lake_only.",
            f"Append inventory summary to {FETCH_LOG_HINT} if desired.",
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
        "policy_risk": list(meta.get("policy_risk") or []),
        "assets": list(meta.get("assets") or []),
        "network_fetch_performed": False,
        "approval_doc": str(APPROVAL_DOC.relative_to(ROOT)),
        "costs_doc": str(COSTS_DOC.relative_to(ROOT)),
        "fetch_log_hint": FETCH_LOG_HINT,
        "manual_recipe": list(meta["manual_recipe"]),
        "fetch_script": "src/v2/extension/scripts/fetch_approved_plan.py",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "note": (
            "Fail-closed stub: no urllib/wget/curl in this CLI. After human approval + env gates, "
            "execute assets via scripts/fetch_approved_plan.py --execute (still requires env match)."
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
            "decision": meta.get("decision"),
            "message": (
                "Download gate closed. Set SCREG_DOWNLOAD_APPROVED=1 and "
                f"SCREG_DOWNLOAD_PLAN_ID={plan_id} after filling "
                f"{APPROVAL_DOC.relative_to(ROOT)}. This CLI never fetches."
                + (
                    " Plan is pending_large (policy-gated infrastructure; not permanently impossible)."
                    if meta.get("decision") == "pending_large"
                    else ""
                )
            ),
            "required_env": {
                APPROVAL_ENV: "1",
                PLAN_MATCH_ENV: plan_id,
            },
            "policy_risk": list(meta.get("policy_risk") or []),
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

    fetch_cmd = (
        f"SCREG_DOWNLOAD_APPROVED=1 SCREG_DOWNLOAD_PLAN_ID={plan_id} "
        f"python src/v2/extension/scripts/fetch_approved_plan.py "
        f"--plan-id {plan_id} --execute"
    )
    plan["next_manual_commands"] = [
        "# Approved dry-run only — CLI performed no network I/O.",
        *[f"#   {step}" for step in plan["manual_recipe"]],
        f"# Real fetch (only after verbal approve + env): {fetch_cmd}",
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
