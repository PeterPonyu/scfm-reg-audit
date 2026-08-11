#!/usr/bin/env python3
"""Emit SI claim-pack tables from existing PeerJ public JSON (zero recompute).

Reads only ``results/*.public.json`` (+ ``paper/panel_data.json`` usability) and
writes under ``docs/reports/extension-claim-pack/`` (MANIFEST-safe local overlay).
Does not touch MANIFEST locks, Support row counts, or the PeerJ submission package.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

EXT_DIR = Path(__file__).resolve().parent
ROOT = EXT_DIR.parents[2]
RESULTS = ROOT / "results"
PANEL_DATA = ROOT / "paper" / "panel_data.json"
OUT_DIR = ROOT / "docs" / "reports" / "extension-claim-pack"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _dual_full(row: dict[str, Any]) -> bool:
    return (
        float(row["mantel"]["bh_q_family"]) < 0.05
        and float(row["degree_preserving"]["bh_q_family"]) < 0.05
    )


def build_dual_null_table(audit: dict[str, Any], oc: dict[str, Any]) -> dict[str, Any]:
    rows = []
    for tissue in ("brain", "pbmc"):
        for row in audit["pooled"][tissue]["rows"]:
            if row.get("row_type") != "pooled_fm" or row.get("confound_spec") != "full":
                continue
            if "mantel" not in row:
                continue
            dual = _dual_full(row)
            rows.append(
                {
                    "tissue": tissue,
                    "model_label": row["model_label"],
                    "rho_full": float(row["observed_partial_rho"]),
                    "q_mantel": float(row["mantel"]["bh_q_family"]),
                    "q_degree": float(row["degree_preserving"]["bh_q_family"]),
                    "dual_null_full": dual,
                    "dual_null_positive_rho": dual and float(row["observed_partial_rho"]) > 0,
                }
            )
    n_dual = sum(1 for r in rows if r["dual_null_full"])
    n_pos = sum(1 for r in rows if r["dual_null_positive_rho"])
    return {
        "schema_version": 1,
        "source": "results/fixed_panel_audit_v2.public.json",
        "oc_source": "results/dual_null_oc_independence_v2.public.json",
        "n_full_rows": len(rows),
        "n_dual_null_full": n_dual,
        "n_dual_null_positive_rho": n_pos,
        "observed_dual_null_full_counts": oc.get("observed_dual_null_full_counts"),
        "independence_oc": {
            tissue: {
                "rate_at_least_one_dual_null_row": fam["rate_at_least_one_dual_null_row"],
                "n_sims": fam["n_sims"],
                "n_rows_in_family": fam["n_rows_in_family"],
            }
            for tissue, fam in (oc.get("families") or {}).items()
        },
        "interpretation": oc.get("interpretation"),
        "rows": rows,
    }


def build_fm_baseline_table(fm: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "source": "results/fm_vs_baseline_observed_v2.public.json",
        "n_fm_rows": fm.get("n_fm_rows"),
        "n_beats_baseline": fm.get("n_beats_baseline"),
        "n_dual_null_full": fm.get("n_dual_null_full"),
        "n_dual_and_beats_baseline": fm.get("n_dual_and_beats_baseline"),
        "summary": fm.get("summary"),
        "rows": fm.get("rows"),
    }


def build_protocol_pass_table(
    audit: dict[str, Any],
    fm: dict[str, Any],
    usability: dict[str, Any],
) -> dict[str, Any]:
    """Reconstruct SAP protocol gates available from public JSON.

    Multi-RO sign is left null (not inventable from a single public field).
    Frozen empirical claim remains 0/13 from SAP / Table 4.
    """
    by_spec: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for tissue in ("brain", "pbmc"):
        for row in audit["pooled"][tissue]["rows"]:
            if row.get("row_type") != "pooled_fm" or "mantel" not in row:
                continue
            by_spec.setdefault((tissue, row["confound_spec"]), {})[row["model_label"]] = row

    fm_index = {(r["tissue"], r["model_label"]): r for r in fm.get("rows") or []}
    out_rows = []
    for (tissue, spec), full_map in list(by_spec.items()):
        if spec != "full":
            continue
        for model_label, full in full_map.items():
            nondeg = by_spec.get((tissue, "non_degree"), {}).get(model_label)
            rho_full = float(full["observed_partial_rho"])
            rho_nd = float(nondeg["observed_partial_rho"]) if nondeg else None
            dual_full = _dual_full(full)
            dual_nd = bool(nondeg) and _dual_full(nondeg)
            conc = usability.get(tissue, {}).get(model_label)
            concordance = None if conc is None else float(conc) > 0
            if rho_nd is None:
                nd_same_sign = None
            else:
                nd_same_sign = (rho_full > 0) == (rho_nd > 0) or (rho_full == 0 and rho_nd == 0)
            beat = fm_index.get((tissue, model_label), {}).get("beats_baseline")
            partial_pass = bool(dual_full and concordance and nd_same_sign and beat)
            out_rows.append(
                {
                    "tissue": tissue,
                    "model_label": model_label,
                    "dual_full": dual_full,
                    "dual_nondeg": dual_nd,
                    "concordance": concordance,
                    "concordance_rho": conc,
                    "nd_same_sign": nd_same_sign,
                    "multi_ro_sign": None,
                    "rho_gt_baseline": beat,
                    "protocol_pass_partial_no_multi_ro": partial_pass,
                    "protocol_pass_frozen_table4": False,
                }
            )

    n_partial = sum(1 for r in out_rows if r["protocol_pass_partial_no_multi_ro"])
    return {
        "schema_version": 1,
        "sap_definition": (
            "Dual full ∧ Concordance ∧ ND same sign ∧ Multi-RO sign ∧ ρ>baseline"
        ),
        "sources": [
            "results/fixed_panel_audit_v2.public.json",
            "results/fm_vs_baseline_observed_v2.public.json",
            "paper/panel_data.json#usability_fm_vs_coexp",
        ],
        "n_full_rows": len(out_rows),
        "n_protocol_pass_frozen": 0,
        "n_protocol_pass_partial_no_multi_ro": n_partial,
        "note": (
            "Frozen PeerJ claim remains 0/13 protocol-pass (SAP §4 / Table 4). "
            "Multi-RO sign is not reconstructed here (null). "
            "Partial conjunction excluding Multi-RO is diagnostic only."
        ),
        "gate_counts": {
            "dual_full": sum(1 for r in out_rows if r["dual_full"]),
            "dual_nondeg": sum(1 for r in out_rows if r["dual_nondeg"]),
            "concordance": sum(1 for r in out_rows if r["concordance"]),
            "nd_same_sign": sum(1 for r in out_rows if r["nd_same_sign"]),
            "rho_gt_baseline": sum(1 for r in out_rows if r["rho_gt_baseline"]),
        },
        "rows": out_rows,
    }


def build_markdown(
    dual: dict[str, Any],
    fm_base: dict[str, Any],
    protocol: dict[str, Any],
) -> str:
    lines = [
        "# Local claim pack (JSON-only, zero download)",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%MZ')}",
        "",
        "Extension SI numbers from existing `results/*.public.json`. "
        "Does **not** change PeerJ Support denominators.",
        "",
        "## Dual-null rarity (C1)",
        "",
        f"- Full-spec FM rows: **{dual['n_full_rows']}**",
        f"- Dual-null full (both BH q < 0.05): **{dual['n_dual_null_full']}**",
        f"- Dual-null with ρ > 0: **{dual['n_dual_null_positive_rho']}**",
        f"- Observed counts by tissue: `{dual['observed_dual_null_full_counts']}`",
        "",
        "Independence OC P(≥1 dual-null row):",
    ]
    for tissue, fam in dual["independence_oc"].items():
        lines.append(
            f"- {tissue}: **{fam['rate_at_least_one_dual_null_row']}** "
            f"(n_sims={fam['n_sims']}, family_n={fam['n_rows_in_family']})"
        )
    lines.extend(
        [
            "",
            f"_{dual.get('interpretation')}_",
            "",
            "## FM − baseline (C5)",
            "",
            f"- {fm_base.get('summary')}",
            f"- Beats baseline: **{fm_base.get('n_beats_baseline')}/{fm_base.get('n_fm_rows')}**",
            f"- Dual ∧ beats baseline: **{fm_base.get('n_dual_and_beats_baseline')}**",
            "",
            "## Protocol-pass gates (C2)",
            "",
            f"- Frozen PeerJ protocol-pass: **{protocol['n_protocol_pass_frozen']}/"
            f"{protocol['n_full_rows']}**",
            f"- Partial conjunction (excl. Multi-RO): "
            f"**{protocol['n_protocol_pass_partial_no_multi_ro']}** "
            "(diagnostic only; not a SAP pass)",
            f"- Gate counts: `{protocol['gate_counts']}`",
            "",
            protocol["note"],
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def emit(out_dir: Path) -> dict[str, Path]:
    audit = load_json(RESULTS / "fixed_panel_audit_v2.public.json")
    oc = load_json(RESULTS / "dual_null_oc_independence_v2.public.json")
    fm = load_json(RESULTS / "fm_vs_baseline_observed_v2.public.json")
    panel = load_json(PANEL_DATA)
    usability = panel["usability_fm_vs_coexp"]

    dual = build_dual_null_table(audit, oc)
    fm_base = build_fm_baseline_table(fm)
    protocol = build_protocol_pass_table(audit, fm, usability)

    if dual["n_dual_null_full"] != fm.get("n_dual_null_full"):
        raise SystemExit(
            f"dual-null count mismatch: audit={dual['n_dual_null_full']} "
            f"fm_vs_baseline={fm.get('n_dual_null_full')}"
        )
    if dual["n_full_rows"] != 13:
        raise SystemExit(f"expected 13 full rows, got {dual['n_full_rows']}")

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "dual_null": out_dir / "dual_null_claim_pack.json",
        "fm_baseline": out_dir / "fm_vs_baseline_claim_pack.json",
        "protocol": out_dir / "protocol_pass_claim_pack.json",
        "markdown": out_dir / "CLAIM_PACK.md",
        "index": out_dir / "index.json",
    }
    paths["dual_null"].write_text(json.dumps(dual, indent=2) + "\n")
    paths["fm_baseline"].write_text(json.dumps(fm_base, indent=2) + "\n")
    paths["protocol"].write_text(json.dumps(protocol, indent=2) + "\n")
    paths["markdown"].write_text(build_markdown(dual, fm_base, protocol))
    index = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "peerj_support_rows_touched": False,
        "files": {k: v.relative_to(ROOT).as_posix() for k, v in paths.items() if k != "index"},
        "counts": {
            "dual_null_full": dual["n_dual_null_full"],
            "beats_baseline": fm_base.get("n_beats_baseline"),
            "protocol_pass_frozen": protocol["n_protocol_pass_frozen"],
        },
        "label_histogram": dict(Counter(r["tissue"] for r in dual["rows"])),
    }
    paths["index"].write_text(json.dumps(index, indent=2) + "\n")
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUT_DIR,
        help="Output directory (default: docs/reports/extension-claim-pack)",
    )
    args = parser.parse_args(argv)
    paths = emit(args.out_dir)
    print("claim pack emitted:")
    for key, path in paths.items():
        print(f"  {key}: {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
