#!/usr/bin/env python3
"""Approval-gated fetch helper for extension download plans.

Unlike ``cli.py download`` / ``download_gate.py`` (always dry-run), this script
*can* perform network I/O — but only when:

1. ``SCREG_DOWNLOAD_APPROVED=1``
2. ``SCREG_DOWNLOAD_PLAN_ID`` matches ``--plan-id``
3. ``--execute`` is passed

Without ``--execute``, it inventories / prints the recipe and exits 0 with
``network_fetch_performed=false``.

D5 defaults to **local inventory only** (RNA lakes already on Desktop); it never
wires assets into Support / ``G_ATAC``.

Do not run D4/D5 ``--execute`` unless a human has verbally approved.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_EXT = Path(__file__).resolve().parents[1]
_ROOT = _EXT.parents[2]
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from download_gate import (  # noqa: E402
    APPROVAL_ENV,
    PLAN_MATCH_ENV,
    PLAN_REGISTRY,
    build_dry_run_plan,
)
from paths import DESKTOP_DATA  # noqa: E402


def _expand(path: str) -> Path:
    return Path(
        os.path.expandvars(path.replace("${DESKTOP_DATA}", str(DESKTOP_DATA)))
    ).expanduser()


def _approval_ok(plan_id: str) -> bool:
    return (
        os.environ.get(APPROVAL_ENV, "").strip() == "1"
        and os.environ.get(PLAN_MATCH_ENV, "").strip() == plan_id
    )


def _inventory_local(dirs: list[str]) -> dict[str, Any]:
    rows = []
    total = 0
    for d in dirs:
        p = _expand(d)
        exists = p.is_dir()
        n_files = 0
        nbytes = 0
        if exists:
            for child in p.iterdir():
                if child.is_file():
                    n_files += 1
                    nbytes += child.stat().st_size
        total += nbytes
        rows.append(
            {
                "path": str(p),
                "exists": exists,
                "n_files": n_files,
                "bytes": nbytes,
            }
        )
    return {"dirs": rows, "total_bytes": total}


def _fetch_url(url: str, dest: Path) -> dict[str, Any]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    # Resume if partial exists.
    existing = tmp.stat().st_size if tmp.exists() else 0
    req = urllib.request.Request(url)
    if existing > 0:
        req.add_header("Range", f"bytes={existing}-")
    with urllib.request.urlopen(req, timeout=120) as resp, open(
        tmp, "ab" if existing else "wb"
    ) as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    return {
        "url": url,
        "dest": str(dest),
        "bytes": dest.stat().st_size,
        "resumed_from": existing,
    }


def run_fetch(
    plan_id: str,
    *,
    execute: bool = False,
    asset_ids: list[str] | None = None,
) -> tuple[int, dict[str, Any]]:
    if plan_id not in PLAN_REGISTRY:
        return 2, {
            "status": "unknown_plan",
            "plan_id": plan_id,
            "known_plan_ids": sorted(PLAN_REGISTRY),
            "network_fetch_performed": False,
        }

    meta = PLAN_REGISTRY[plan_id]
    if meta.get("rejected"):
        return 2, {
            "status": "rejected_plan",
            "plan_id": plan_id,
            "message": f"Plan {plan_id} is permanently rejected.",
            "network_fetch_performed": False,
        }

    if not _approval_ok(plan_id):
        return 2, {
            "status": "approval_required",
            "plan_id": plan_id,
            "message": (
                "Fetch script closed. Set SCREG_DOWNLOAD_APPROVED=1 and "
                f"SCREG_DOWNLOAD_PLAN_ID={plan_id}. Without --execute this script "
                "still only inventories / dry-runs."
            ),
            "required_env": {APPROVAL_ENV: "1", PLAN_MATCH_ENV: plan_id},
            "network_fetch_performed": False,
        }

    plan = build_dry_run_plan(plan_id)
    assets = list(meta.get("assets") or [])
    if asset_ids:
        want = set(asset_ids)
        assets = [a for a in assets if a.get("id") in want]

    if not execute:
        inventories = []
        for asset in assets:
            if asset.get("fetch_via") == "local_inventory" or asset.get("local_dirs"):
                inventories.append(
                    {
                        "id": asset.get("id"),
                        "inventory": _inventory_local(list(asset.get("local_dirs") or [])),
                    }
                )
        payload = {
            **plan,
            "status": "approved_fetch_dry_run",
            "execute": False,
            "selected_assets": [a.get("id") for a in assets],
            "local_inventory": inventories,
            "network_fetch_performed": False,
            "message": (
                "Dry-run only. Re-run with --execute after verbal approval to fetch "
                "HTTP assets (D4) or re-inventory (D5)."
            ),
            "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        }
        return 0, payload

    # Execute path.
    results: list[dict[str, Any]] = []
    network = False
    for asset in assets:
        aid = asset.get("id")
        via = asset.get("fetch_via")
        if via == "local_inventory" or asset.get("local_dirs"):
            results.append(
                {
                    "id": aid,
                    "action": "local_inventory",
                    "inventory": _inventory_local(list(asset.get("local_dirs") or [])),
                    "network_fetch_performed": False,
                    "g_atac_forbidden": plan_id == "D5",
                }
            )
            continue
        if via == "synapse_manual" or not asset.get("url"):
            results.append(
                {
                    "id": aid,
                    "action": "skipped_manual",
                    "reason": asset.get("checksum_note")
                    or "No auto-URL; fetch manually (Synapse / human shortlist).",
                    "network_fetch_performed": False,
                }
            )
            continue
        dest = _expand(str(asset["dest"]))
        fetched = _fetch_url(str(asset["url"]), dest)
        fetched["id"] = aid
        fetched["action"] = "fetched"
        fetched["network_fetch_performed"] = True
        results.append(fetched)
        network = True

    payload = {
        **plan,
        "status": "executed",
        "execute": True,
        "results": results,
        "network_fetch_performed": network,
        "g_atac_forbidden_reminder": (
            "D5 RNA lakes must not enter Support / G_ATAC even after inventory/fetch."
            if plan_id == "D5"
            else "D4 lakes stay extension/lake_blocked until subsetting; never PeerJ Support."
            if plan_id == "D4"
            else None
        ),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
    }
    return 0, payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Approval-gated plan fetch. Requires SCREG_DOWNLOAD_APPROVED=1 + "
            "matching SCREG_DOWNLOAD_PLAN_ID. Default is dry-run; pass --execute "
            "only after verbal human approval."
        )
    )
    parser.add_argument("--plan-id", required=True, help="D0–D5 plan id")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Perform HTTP fetch / local inventory (still requires env approval)",
    )
    parser.add_argument(
        "--asset-id",
        action="append",
        default=[],
        help="Optional asset id filter (repeatable); default = all plan assets",
    )
    args = parser.parse_args(argv)
    code, payload = run_fetch(
        args.plan_id,
        execute=args.execute,
        asset_ids=args.asset_id or None,
    )
    print(json.dumps(payload, indent=2))
    if code != 0:
        print(payload.get("message", "fetch refused"), file=sys.stderr)
    elif not args.execute:
        print(
            "Exit 0 = dry-run / inventory recipe only — no download claimed.",
            file=sys.stderr,
        )
    else:
        print(
            "Execute finished (see JSON). network_fetch_performed="
            f"{payload.get('network_fetch_performed')}",
            file=sys.stderr,
        )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
