#!/usr/bin/env python3
"""Construct-lane stubs (fibroblast-style Mantel / decomp hooks).

These helpers document the env seam used by ``build_atac_graph_v2.py`` and
emit extension-only dry-run plans. They do **not** build new Support rows or
write PeerJ public JSON.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))

from registry import HEAVY_ARTIFACT_ROOT, load_extension_registry  # noqa: E402

ROOT = _EXT.parents[2]
PLAN_ROOT = ROOT / "docs" / "reports" / "extension-plans" / "construct"


def fibro_style_env(tissue_id: str, atac_file: str | None = None) -> dict[str, str]:
    """Return env vars mirroring the GSE206767 construct seam."""
    reg = load_extension_registry()
    reg.assert_may_emit_g_atac(tissue_id)
    meta = reg.get_tissue(tissue_id)
    tag = str(meta.get("g_atac_tag") or tissue_id)
    # Keep paths relative in emitted plans so capsule privacy checks stay green.
    env = {
        "TAG": tag,
        "ATAC_FILE": atac_file or os.environ.get("ATAC_FILE", ""),
        "SCREG_EXTENSION_OUT": f"results/v2/extension/construct/{tag}",
        "SCREG_EXTENSION_LANE": "construct",
        "SCREG_PEERJ_SUPPORT_LOCK": "1",
    }
    return env


def mantel_decomp_plan(tissue_id: str) -> dict[str, Any]:
    """Plan Mantel / additive-decomp vs locked brain+PBMC+fibro proxies."""
    reg = load_extension_registry()
    reg.assert_may_emit_g_atac(tissue_id)
    meta = reg.get_tissue(tissue_id)
    tag = meta.get("g_atac_tag") or tissue_id
    return {
        "tissue_id": tissue_id,
        "lane": "construct",
        "g_atac_tag": tag,
        "compare_to_tags": ["GSE174367", "PBMC10k", "GSE206767"],
        "outputs": {
            "g_atac": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/G_ATAC_v2_{tag}.npz",
            "mantel": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/mantel_vs_locked.json",
            "decomp": f"{HEAVY_ARTIFACT_ROOT}construct/{tag}/additive_decomp_row.json",
        },
        "peerj_decomp_rows_unchanged": 3,
        "peerj_support_rows_unchanged": 13,
        "status": "stub",
        "next_steps": [
            "Provide ATAC_FILE + cell-type META_FILE for the tissue",
            "Run build_atac_graph_v2.py with TAG set; redirect NPZ under results/v2/extension/",
            "Emit Mantel/decomp into results/v2/extension/ only (do not rewrite "
            "cross_tissue_additive_decomp_v2.public.json)",
        ],
        "panel_policy": meta.get("panel_policy"),
        "notes": meta.get("notes"),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissue", default="bmmc")
    parser.add_argument("--atac-file", default="")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write dry-run JSON under docs/reports/extension-plans/",
    )
    args = parser.parse_args(argv)

    plan = mantel_decomp_plan(args.tissue)
    plan["env"] = fibro_style_env(args.tissue, atac_file=args.atac_file or None)
    print(json.dumps(plan, indent=2))

    if args.write:
        out = PLAN_ROOT / str(plan["g_atac_tag"])
        out.mkdir(parents=True, exist_ok=True)
        path = out / "construct_plan.dry_run.json"
        path.write_text(json.dumps(plan, indent=2) + "\n")
        print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
