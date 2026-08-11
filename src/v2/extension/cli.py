#!/usr/bin/env python3
"""First-class CLI for extension-lane infrastructure (local assets only)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_EXT = Path(__file__).resolve().parent
if str(_EXT) not in sys.path:
    sys.path.insert(0, str(_EXT))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="extension",
        description=(
            "scReg-Eval extension infra "
            "(registry / claim-pack / construct / baselines / fail-closed download)"
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_reg = sub.add_parser("registry", help="Load tissue/method registry + deny-lake checks")
    p_reg.add_argument("--dry-run-tissue", default="bmmc")
    p_reg.add_argument("--json", action="store_true")

    p_claim = sub.add_parser("claim-pack", help="Emit SI claim-pack tables from public JSON")
    p_claim.add_argument("--out-dir", type=Path, default=None)

    p_con = sub.add_parser("construct", help="Construct-lane Mantel/decomp hooks")
    p_con.add_argument("--tissue", default="fibroblast")
    p_con.add_argument("--atac-file", default="")
    p_con.add_argument("--execute", action="store_true")
    p_con.add_argument("--write", action="store_true")

    p_base = sub.add_parser("baselines", help="Emit Tier A–C baseline artifacts")
    p_base.add_argument("--method", default="motif_only_rp")
    p_base.add_argument("--all", action="store_true")
    p_base.add_argument("--execute", action="store_true")
    p_base.add_argument("--proxy-tag", default="GSE174367")
    p_base.add_argument("--write", action="store_true")

    p_dl = sub.add_parser(
        "download",
        help=(
            "Fail-closed download gate (no network). "
            "Exit 0 = recipe/dry-run only, never download success. "
            "SCREG_DOWNLOAD_* env is checklist ceremony, not auth."
        ),
    )
    p_dl.add_argument("--plan-id", required=True, help="Approval matrix id (D0–D5)")
    p_dl.add_argument(
        "--no-write",
        action="store_true",
        help="Skip writing dry-run JSON under results/v2/extension/download-plans/",
    )

    p_d1 = sub.add_parser(
        "descartes-bridge",
        help="D1 local RDS/h5ad → ATAC_FILE readiness (no fetch)",
    )
    p_d1.add_argument("--pilot-dir", default="")
    p_d1.add_argument("--write", action="store_true")

    p_bmmc = sub.add_parser(
        "prepare-bmmc",
        help="Extract BMMC ATAC peaks into extension overlay (no fetch / no G_ATAC)",
    )
    p_bmmc.add_argument("--src", default="")
    p_bmmc.add_argument("--execute", action="store_true")
    p_bmmc.add_argument("--write", action="store_true")

    args = parser.parse_args(argv)

    if args.cmd == "registry":
        import registry as reg

        return reg.main(
            ["--dry-run-tissue", args.dry_run_tissue] + (["--json"] if args.json else [])
        )
    if args.cmd == "claim-pack":
        import emit_claim_pack as ecp

        argv2: list[str] = []
        if args.out_dir is not None:
            argv2.extend(["--out-dir", str(args.out_dir)])
        return ecp.main(argv2)
    if args.cmd == "construct":
        import construct_hooks as ch

        argv2 = ["--tissue", args.tissue]
        if args.atac_file:
            argv2.extend(["--atac-file", args.atac_file])
        if args.execute:
            argv2.append("--execute")
        if args.write:
            argv2.append("--write")
        return ch.main(argv2)
    if args.cmd == "baselines":
        import baseline_stubs as bs

        argv2 = ["--method", args.method, "--proxy-tag", args.proxy_tag]
        if args.all:
            argv2.append("--all")
        if args.execute:
            argv2.append("--execute")
        if args.write:
            argv2.append("--write")
        return bs.main(argv2)
    if args.cmd == "download":
        import download_gate as dg

        argv2 = ["--plan-id", args.plan_id]
        if args.no_write:
            argv2.append("--no-write")
        return dg.main(argv2)
    if args.cmd == "descartes-bridge":
        import descartes_bridge as db

        argv2: list[str] = []
        if args.pilot_dir:
            argv2.extend(["--pilot-dir", args.pilot_dir])
        if args.write:
            argv2.append("--write")
        return db.main(argv2)
    if args.cmd == "prepare-bmmc":
        import bmmc_prepare as bp

        argv2 = []
        if args.src:
            argv2.extend(["--src", args.src])
        if args.execute:
            argv2.append("--execute")
        if args.write:
            argv2.append("--write")
        return bp.main(argv2)
    raise AssertionError(f"unhandled cmd: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
