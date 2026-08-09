#!/usr/bin/env python3
"""Validate the sanitized scReg-Eval audit capsule artifacts (v0.2.0 contract)."""
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"

# Order must match first-to-last \\input{figs/...} in paper/manuscript.tex.
CURRENT_FIGURES = (
    "fig10_coverage_qc.tex",
    "fig1_truth_construct.tex",
    "fig2_cross_tissue_decomp.tex",
    "fig3_primary_audit.tex",
    "fig4_usability_check.tex",
    "fig5_null_diagnostics.tex",
    "fig6_spec_sensitivity.tex",
    "fig7_pertype_descriptive.tex",
    "fig8_injection_ladder.tex",
    "fig9_tf_probe.tex",
    "fig11_third_tissue_transfer.tex",
    "fig12_protocol_pass_matrix.tex",
    "fig13_scope_card.tex",
)
CURRENT_TABLES = (
    "table5_related_work.tex",
    "table2_cross_tissue_observed.tex",
    "table1_primary_fixed_panel.tex",
    "table4_protocol_pass.tex",
    "table6_fm_vs_baseline.tex",
    "table3_pertype_ranges.tex",
)
CURRENT_FRAGMENTS = CURRENT_FIGURES + CURRENT_TABLES
FIGURE_INPUT_RE = re.compile(r"\\input\{figs/([^}]+\.tex)\}")


def load(name):
    return json.loads((RESULTS / name).read_text(),
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def walk_numbers(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from walk_numbers(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_numbers(child)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield value


def bh(pvalues):
    n = len(pvalues)
    order = sorted(range(n), key=pvalues.__getitem__)
    out = [0.0] * n
    running = 1.0
    for reverse_index in range(n - 1, -1, -1):
        index = order[reverse_index]
        running = min(running, pvalues[index] * n / (reverse_index + 1), 1.0)
        out[index] = running
    return out


def check_no_private_paths():
    # The validator necessarily names the forbidden needles; it must not flag itself.
    forbidden = ["/home/zeyufu", "zeyufu/Desktop", ".omc/", "__pycache__"]
    for path in ROOT.rglob("*"):
        if path.name == "validate_artifacts.py":
            continue
        if path.is_file() and path.suffix in {".json", ".md", ".tex", ".py", ".R", ".bib", ".cff", ".txt"}:
            text = path.read_text(errors="replace")
            for needle in forbidden:
                assert needle not in text, f"private path {needle} in {path.relative_to(ROOT)}"


def figure_inputs(path):
    return tuple(FIGURE_INPUT_RE.findall(path.read_text()))


def check_figure_contract():
    manuscript_inputs = figure_inputs(ROOT / "paper/manuscript.tex")
    manuscript_figures = tuple(name for name in manuscript_inputs if name.startswith("fig"))
    manuscript_tables = tuple(name for name in manuscript_inputs if name.startswith("table"))
    preview_inputs = figure_inputs(ROOT / "paper/figs_preview.tex")
    bundled = {path.name for path in (ROOT / "paper/figs").glob("*.tex")}

    assert manuscript_figures == CURRENT_FIGURES
    assert manuscript_tables == CURRENT_TABLES
    assert preview_inputs == CURRENT_FRAGMENTS
    assert bundled == set(CURRENT_FRAGMENTS), (
        f"active figure allowlist mismatch: {sorted(bundled)}")


def main():
    audit = load("fixed_panel_audit_v2.public.json")
    injection = load("fixed_panel_signal_injection_v2.public.json")
    status = load("inference_status_v2.public.json")
    brain_base = load("brain_coexp_baseline_null_v2.public.json")
    pbmc_base = load("pbmc_coexp_baseline_null_v2.public.json")
    probe = load("tf_probe_pair_stats_v2.public.json")
    sensitivity = load("tf_probe_pair_sensitivity_v2.public.json")
    invariance = load("proxy_celltype_invariance_v2.public.json")
    subdivided = load("injection_subdivided_v2.public.json")
    decomp = load("cross_tissue_additive_decomp_v2.public.json")

    docs = [audit, injection, status, brain_base, pbmc_base, probe, sensitivity, invariance,
            subdivided, decomp]
    assert all(math.isfinite(float(x)) for doc in docs for x in walk_numbers(doc))
    assert audit["panel"]["n_tf"] == 446 and audit["panel"]["n_genes"] == 1200
    assert len(decomp["rows"]) == 3

    rows = []
    for tissue in ("brain", "pbmc"):
        rows.extend(row for row in audit["pooled"][tissue]["rows"] if row.get("row_type") == "pooled_fm")
    assert len(rows) == 26, f"expected 26 pooled rows (13 per spec), got {len(rows)}"
    full_rows = [row for row in rows if row["confound_spec"] == "full"]
    assert len(full_rows) == 13, f"expected 13 full-spec rows, got {len(full_rows)}"

    families = defaultdict(list)
    for row in rows:
        for key in ("mantel", "degree_preserving"):
            result = row[key]
            assert result["N_perm"] == 999 and result["resolution"] == 0.001
            expected = (result["null_obs_count_at_or_above_obs"] + 1) / 1000
            assert abs(result["p_mc"] - expected) < 1e-12
            families[result["family_id"]].append(result)
    assert len(families) == 8 and set(families) == set(status["BH_family_definitions"])
    for results in families.values():
        expected = bh([result["p_mc"] for result in results])
        assert all(abs(result["bh_q_family"] - q) < 1e-6 for result, q in zip(results, expected))

    pertype = []
    for tissue in ("brain", "pbmc"):
        for spec in ("full_confound", "non_degree_confound"):
            pertype.extend(row for row in audit["per_cell_type"][tissue][spec]["rows"]
                           if row.get("row_type") == "pertype_fm")
    assert len(pertype) == 58, f"expected 58 per-type rows, got {len(pertype)}"
    assert all(not any(key in row for key in ("p_mc", "bh_q_family", "mantel", "degree_preserving"))
               for row in pertype)
    assert len(audit["cross_tissue_construct_reproducibility"]["rows"]) == 3

    alphas = [0.0, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.75, 1.0]
    seeds = []
    for tissue in ("brain", "pbmc"):
        tissue_rows = injection["tissues"][tissue]["rows"]
        assert [row["alpha"] for row in tissue_rows] == alphas
        for row in tissue_rows:
            assert len(row["replicate_runs"]) == 30
            seeds.extend(run["seed"] for run in row["replicate_runs"])
    assert len(seeds) == 660 and len(set(seeds)) == 660

    assert brain_base["seed_contract"] == "explicit_integer_v1"
    assert brain_base["n_perm"] == 999
    assert probe["schema_version"] == 2
    expected_families = {"co_expression", "geneformer_embed", "geneformer_attn",
                         "scGPT_encoder", "UCE_encoder", "random_floor"}
    assert set(probe["families"]) == expected_families
    assert set(probe["contrasts_vs_baseline"]) == expected_families - {"co_expression"}
    assert all("mantel_seed" in fam for fam in probe["families"].values())
    assert probe["n_perm"] == 999

    manifest = json.loads((ROOT / "MANIFEST.json").read_text())
    records = manifest["files"]
    listed = {record["path"] for record in records}
    assert len(listed) == len(records)
    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.name != "MANIFEST.json"
        and "__pycache__" not in path.parts and ".omc" not in path.parts
        and ".pytest_cache" not in path.parts and path.suffix != ".pyc"
    }
    assert listed == actual, (f"manifest coverage mismatch: "
                              f"missing={sorted(actual - listed)}, extra={sorted(listed - actual)}")
    for record in records:
        path = ROOT / record["path"]
        assert record["bytes"] == path.stat().st_size
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]

    check_figure_contract()
    check_no_private_paths()
    version = manifest.get("version", "unknown")
    print(f"PASS: scReg-Eval audit capsule {version} artifacts, figures, manifest, and privacy are consistent")


if __name__ == "__main__":
    sys.exit(main())
