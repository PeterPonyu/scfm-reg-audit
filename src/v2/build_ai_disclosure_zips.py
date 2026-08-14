#!/usr/bin/env python3
"""Build the three PeerJ CS AI-in-code supplemental zip archives.

PeerJ Author Instructions name supplemental files as
'Supplemental [Item] S[number]'. Archives are written under
paper/submission_peerj/supplemental/, never into flat_upload/, so they cannot
leak into the manuscript checksum set or GitHub Pages assemble output.
"""
import hashlib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "paper" / "submission_peerj"
SRC = PKG / "ai_code_disclosure"
OUT = PKG / "supplemental"

ARCHIVES = (
    ("Supplemental_Data_S1_original_code.zip",
     "Supplemental Data S1",
     "Original (pre-engineering) copies of the three analysis scripts.",
     ("before/pbmc_uce_eval_v2.py",
      "before/pair_probe_stats.py",
      "before/fixed_panel_audit.py")),
    ("Supplemental_Data_S2_prompts.zip",
     "Supplemental Data S2",
     "Reconstructed prompt log for the three disclosed engineering edits.",
     ("prompts/00_INDEX.md",)),
    ("Supplemental_Data_S3_resulting_code.zip",
     "Supplemental Data S3",
     "Resulting (post-engineering) copies of the three analysis scripts.",
     ("after/pbmc_uce_eval_v2.py",
      "after/pair_probe_stats.py",
      "after/fixed_panel_audit.py")),
)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    if not SRC.is_dir():
        raise RuntimeError(f"missing disclosure tree: {SRC}")
    OUT.mkdir(parents=True, exist_ok=True)
    for stale in PKG.glob("ai_*.zip"):
        stale.unlink()
    checksum_lines = []
    for zip_name, si_label, purpose, members in ARCHIVES:
        zip_path = OUT / zip_name
        if zip_path.exists():
            zip_path.unlink()
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            readme = (
                "scReg-Eval PeerJ CS AI-in-code supplemental archive\n"
                f"{si_label}\n"
                f"File: {zip_name}\n"
                f"{purpose}\n"
                "Auxiliary engineering disclosure only.\n"
                "Companion archives: Supplemental Data S1 (original code), "
                "S2 (prompts), S3 (resulting code).\n"
                "Production code remains src/v2/ in the repository.\n"
            )
            zf.writestr("README.txt", readme)
            for rel in members:
                src = SRC / rel
                if not src.is_file():
                    raise RuntimeError(f"missing disclosure member: {rel}")
                zf.write(src, arcname=rel)
        checksum_lines.append(f"{sha256_file(zip_path)}  {zip_name}")
        print(f"wrote {zip_path} ({zip_path.stat().st_size} bytes)")
    sums = OUT / "SHA256SUMS_AI.txt"
    sums.write_text("\n".join(checksum_lines) + "\n")
    pkg_sums = PKG / "SHA256SUMS_AI.txt"
    if pkg_sums.exists():
        pkg_sums.unlink()
    print(f"wrote {sums}")


if __name__ == "__main__":
    build()
