#!/usr/bin/env python3
"""Load and validate extension tissue/method registries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

EXT_DIR = Path(__file__).resolve().parent
ROOT = EXT_DIR.parents[2]
DEFAULT_TISSUES = EXT_DIR / "configs" / "tissues.json"
DEFAULT_METHODS = EXT_DIR / "configs" / "methods.json"
FORBIDDEN_G_ATAC_ROLES = frozenset(
    {"out_of_scope", "rna_lake_only", "lake_blocked"}
)
# Heavy extension artifacts under results/v2/extension/ (PeerJ freeze untouched).
# SI claim-pack tables use docs/reports/extension-claim-pack/ so MANIFEST stays frozen.
try:
    from paths import CLAIM_PACK_ROOT, HEAVY_ARTIFACT_ROOT
except ImportError:  # pragma: no cover - script/module dual import
    HEAVY_ARTIFACT_ROOT = "results/v2/extension/"
    CLAIM_PACK_ROOT = "docs/reports/extension-claim-pack/"


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object at top level")
    return data


class ExtensionRegistry:
    """Config-driven tissue + method registry with PeerJ/extension fences."""

    def __init__(self, tissues_doc: dict[str, Any], methods_doc: dict[str, Any]):
        self.tissues_doc = tissues_doc
        self.methods_doc = methods_doc
        self.tissues: dict[str, dict[str, Any]] = dict(tissues_doc.get("tissues") or {})
        self.methods: dict[str, dict[str, Any]] = dict(methods_doc.get("methods") or {})

    def get_tissue(self, tissue_id: str) -> dict[str, Any]:
        try:
            return self.tissues[tissue_id]
        except KeyError as exc:
            raise KeyError(f"unknown tissue id: {tissue_id}") from exc

    def get_method(self, method_id: str) -> dict[str, Any]:
        try:
            return self.methods[method_id]
        except KeyError as exc:
            raise KeyError(f"unknown method id: {method_id}") from exc

    def peerj_tissues(self) -> list[str]:
        return sorted(tid for tid, meta in self.tissues.items() if meta.get("peerj_freeze"))

    def extension_tissues(self, lane: str | None = None) -> list[str]:
        out = []
        for tid, meta in self.tissues.items():
            if meta.get("peerj_freeze"):
                continue
            if lane is not None and meta.get("lane") != lane:
                continue
            out.append(tid)
        return sorted(out)

    def assert_may_emit_g_atac(self, tissue_id: str) -> None:
        meta = self.get_tissue(tissue_id)
        role = meta.get("role")
        if role in FORBIDDEN_G_ATAC_ROLES or meta.get("allow_g_atac") is False:
            raise PermissionError(
                f"tissue {tissue_id!r} role={role!r} is forbidden for G_ATAC / Support"
            )
        if meta.get("lane") not in {"construct", "audit"}:
            raise PermissionError(
                f"tissue {tissue_id!r} lane={meta.get('lane')!r} cannot emit G_ATAC"
            )

    def dry_run_register(self, tissue_id: str) -> dict[str, Any]:
        """Return a dry-run payload for an extension tissue (no PeerJ writes)."""
        meta = self.get_tissue(tissue_id)
        self.assert_may_emit_g_atac(tissue_id)
        return {
            "tissue_id": tissue_id,
            "lane": meta.get("lane"),
            "role": meta.get("role"),
            "panel_policy": meta.get("panel_policy"),
            "peerj_freeze": bool(meta.get("peerj_freeze")),
            "artifact_root": HEAVY_ARTIFACT_ROOT,
            "claim_pack_root": CLAIM_PACK_ROOT,
            "peerj_support_rows_touched": False,
            "notes": meta.get("notes"),
        }

    def summary(self) -> dict[str, Any]:
        return {
            "n_tissues": len(self.tissues),
            "n_methods": len(self.methods),
            "peerj_tissues": self.peerj_tissues(),
            "extension_construct": self.extension_tissues(lane="construct"),
            "extension_audit_candidates": [
                tid
                for tid in self.extension_tissues()
                if self.tissues[tid].get("role") == "primary_audit"
            ],
            # construct_candidate (e.g. BMMC) is intentionally excluded until re-SAP
            "extension_construct_candidates": [
                tid
                for tid in self.extension_tissues()
                if self.tissues[tid].get("role")
                in {"construct", "construct_candidate"}
            ],
            "forbidden_g_atac": sorted(
                tid
                for tid, meta in self.tissues.items()
                if meta.get("role") in FORBIDDEN_G_ATAC_ROLES
                or meta.get("allow_g_atac") is False
            ),
        }


def load_extension_registry(
    tissues_path: Path | None = None,
    methods_path: Path | None = None,
) -> ExtensionRegistry:
    tissues = _load_json(tissues_path or DEFAULT_TISSUES)
    methods = _load_json(methods_path or DEFAULT_METHODS)
    return ExtensionRegistry(tissues, methods)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run-tissue", default="bmmc", help="Tissue id for dry-run")
    parser.add_argument("--json", action="store_true", help="Emit JSON summary")
    args = parser.parse_args(argv)
    reg = load_extension_registry()
    payload = {
        "summary": reg.summary(),
        "dry_run": reg.dry_run_register(args.dry_run_tissue),
    }
    blocked = []
    for tid in (
        "cancer_rna_lakes",
        "development_rna_lakes",
        "descartes_whole_lake",
        "htan_whole_lake",
    ):
        try:
            reg.assert_may_emit_g_atac(tid)
            blocked.append({"tissue_id": tid, "policy": "UNEXPECTED_ALLOW"})
        except PermissionError as exc:
            blocked.append({"tissue_id": tid, "policy": "blocked", "error": str(exc)})
    payload["negative_g_atac_policy"] = blocked
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("extension registry OK")
        print(f"  peerj tissues: {payload['summary']['peerj_tissues']}")
        print(f"  construct candidates: {payload['summary']['extension_construct']}")
        print(f"  dry-run {args.dry_run_tissue}: lane={payload['dry_run']['lane']}")
        print(f"  RNA-lake G_ATAC policy: {[b['policy'] for b in blocked]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
