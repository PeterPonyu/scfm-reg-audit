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
EXCLUDED_PARTS = {"__pycache__", ".omc", ".omx", ".pytest_cache", ".git", "archive"}
EXCLUDED_NAMES = {"MANIFEST.json", ".gitignore"}
# Local/dev trees may coexist with the slim public capsule. These prefixes are
# ignored by the closed-tree equality check; MANIFEST records are still
# hash-verified unconditionally.
# Paths allowed on a developer checkout of the capsule without being listed in
# MANIFEST.json (local overlays, PeerJ rebuild tree, heavy results, agent state).
LOCAL_WORKTREE_PREFIXES = (
    "src/v2/",
    "paper/submission_peerj/",
    "paper/docs/",
    "results/v2/",
    "paper/.tikz",
    ".cursor/",
    ".grok/",
    ".omc/",
    ".omx/",
    "docs/reports/",
    "archive/",
)
LOCAL_WORKTREE_NAMES = {
    "PAPER_REVIEW_TARGETS.md",
    "make_figures.py",
    "wlpeerj.cls",
    "manuscript.bbl",
    "manuscript.docx",
    ".tikz_metrics_pdftex",
    # LaTeX intermediates (never capsule content)
    "manuscript.aux",
    "manuscript.log",
    "manuscript.out",
    "manuscript.fls",
    "manuscript.fdb_latexmk",
    "manuscript.blg",
    "Rplots.pdf",
}

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
    "table3_pertype_ranges.tex",
)
CURRENT_FRAGMENTS = CURRENT_FIGURES + CURRENT_TABLES
FIGURE_INPUT_RE = re.compile(r"\\input\{figs/([^}]+\.tex)\}")


class ValidationError(Exception):
    """A capsule artifact violated the release contract."""


def require(condition, message):
    """Contract gate that survives `python -O`, unlike a bare assert."""
    if not condition:
        raise ValidationError(message)


def load(name):
    def reject_constant(value):
        raise ValidationError(f"{name} contains the non-JSON constant {value}")

    return json.loads((RESULTS / name).read_text(), parse_constant=reject_constant)


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
    # The validator necessarily names the forbidden needles; it must not flag itself
    # or ENVIRONMENT.example, which documents those needles for reviewers.
    forbidden = ["/home/zeyufu", "zeyufu/Desktop"]  # directory caches (.omc/__pycache__) excluded via EXCLUDED_PARTS
    allow_needle_docs = {
        "validate_artifacts.py",
        "ENVIRONMENT.example",
        "validate_capsule.py",
        "test_validate_artifacts.py",  # constructs synthetic leak fixtures
        "build_release_capsule.py",  # documents PATH_REPLACEMENTS scrub needles
    }
    for path in ROOT.rglob("*"):
        if path.name in allow_needle_docs:
            continue
        if EXCLUDED_PARTS.intersection(path.parts):
            continue
        rel = path.relative_to(ROOT).as_posix()
        # Full private result trees may contain machine paths; they are outside the
        # public capsule and are ignored by closed-tree equality. Skip scrub there.
        if rel.startswith("results/v2/"):
            continue
        if path.is_file() and path.suffix in {".json", ".md", ".tex", ".py", ".R", ".bib", ".cff", ".txt"}:
            text = path.read_text(errors="replace")
            for needle in forbidden:
                require(needle not in text,
                        f"private path {needle} in {rel}")


def figure_inputs(path):
    return tuple(FIGURE_INPUT_RE.findall(path.read_text()))


def check_figure_contract():
    manuscript_inputs = figure_inputs(ROOT / "paper/manuscript.tex")
    manuscript_figures = tuple(name for name in manuscript_inputs if name.startswith("fig"))
    manuscript_tables = tuple(name for name in manuscript_inputs if name.startswith("table"))
    preview_inputs = figure_inputs(ROOT / "paper/figs_preview.tex")
    bundled = {path.name for path in (ROOT / "paper/figs").glob("*.tex")}

    require(manuscript_figures == CURRENT_FIGURES,
            f"manuscript figure inputs mismatch: {list(manuscript_figures)}")
    require(manuscript_tables == CURRENT_TABLES,
            f"manuscript table inputs mismatch: {list(manuscript_tables)}")
    require(preview_inputs == CURRENT_FRAGMENTS,
            f"figs_preview.tex inputs mismatch: {list(preview_inputs)}")
    require(bundled == set(CURRENT_FRAGMENTS),
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
    require(all(math.isfinite(float(x)) for doc in docs for x in walk_numbers(doc)),
            "an artifact JSON contains a non-finite number")
    require(audit["panel"]["n_tf"] == 446 and audit["panel"]["n_genes"] == 1200,
            f"panel must be 446 TF x 1200 genes; got {audit['panel']}")
    require(len(decomp["rows"]) == 3,
            f"expected 3 cross-tissue decomposition rows, got {len(decomp['rows'])}")

    rows = []
    for tissue in ("brain", "pbmc"):
        rows.extend(row for row in audit["pooled"][tissue]["rows"] if row.get("row_type") == "pooled_fm")
    require(len(rows) == 26, f"expected 26 pooled rows (13 per spec), got {len(rows)}")
    full_rows = [row for row in rows if row["confound_spec"] == "full"]
    require(len(full_rows) == 13, f"expected 13 full-spec rows, got {len(full_rows)}")

    families = defaultdict(list)
    for row in rows:
        for key in ("mantel", "degree_preserving"):
            result = row[key]
            require(result["N_perm"] == 999 and result["resolution"] == 0.001,
                    f"{row['model_label']}:{key} must use 999 perms at resolution 0.001")
            expected = (result["null_obs_count_at_or_above_obs"] + 1) / 1000
            require(abs(result["p_mc"] - expected) < 1e-12,
                    f"{row['model_label']}:{key} p_mc {result['p_mc']} != plus-one value {expected}")
            families[result["family_id"]].append(result)
    require(len(families) == 8 and set(families) == set(status["BH_family_definitions"]),
            f"BH families mismatch: {sorted(families)}")
    for family_id, results in families.items():
        expected = bh([result["p_mc"] for result in results])
        require(all(abs(result["bh_q_family"] - q) < 1e-6 for result, q in zip(results, expected)),
                f"BH q-values in family {family_id} disagree with recomputation")

    pertype = []
    for tissue in ("brain", "pbmc"):
        for spec in ("full_confound", "non_degree_confound"):
            pertype.extend(row for row in audit["per_cell_type"][tissue][spec]["rows"]
                           if row.get("row_type") == "pertype_fm")
    require(len(pertype) == 58, f"expected 58 per-type rows, got {len(pertype)}")
    require(all(not any(key in row for key in ("p_mc", "bh_q_family", "mantel", "degree_preserving"))
                for row in pertype),
            "per-type rows are descriptive only and must carry no inferential keys")
    require(len(audit["cross_tissue_construct_reproducibility"]["rows"]) == 3,
            "expected 3 cross-tissue construct reproducibility rows")

    alphas = [0.0, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.75, 1.0]
    seeds = []
    for tissue in ("brain", "pbmc"):
        tissue_rows = injection["tissues"][tissue]["rows"]
        require([row["alpha"] for row in tissue_rows] == alphas,
                f"{tissue} injection alphas mismatch: {[row['alpha'] for row in tissue_rows]}")
        for row in tissue_rows:
            require(len(row["replicate_runs"]) == 30,
                    f"{tissue} alpha={row['alpha']} has {len(row['replicate_runs'])} replicates, expected 30")
            seeds.extend(run["seed"] for run in row["replicate_runs"])
    require(len(seeds) == 660 and len(set(seeds)) == 660,
            f"injection seeds must be 660 distinct values; got {len(seeds)} ({len(set(seeds))} distinct)")

    require(brain_base["seed_contract"] == "explicit_integer_v1",
            f"brain baseline seed contract is {brain_base['seed_contract']!r}")
    require(brain_base["n_perm"] == 999, f"brain baseline n_perm is {brain_base['n_perm']}")
    require(probe["schema_version"] == 2, f"probe schema_version is {probe['schema_version']}")
    expected_families = {"co_expression", "geneformer_embed", "geneformer_attn",
                         "scGPT_encoder", "UCE_encoder", "random_floor"}
    require(set(probe["families"]) == expected_families,
            f"probe families mismatch: {sorted(probe['families'])}")
    require(set(probe["contrasts_vs_baseline"]) == expected_families - {"co_expression"},
            f"probe contrasts mismatch: {sorted(probe['contrasts_vs_baseline'])}")
    require(all("mantel_seed" in fam for fam in probe["families"].values()),
            "every probe family must record its mantel_seed")
    require(probe["n_perm"] == 999, f"probe n_perm is {probe['n_perm']}")

    manifest = json.loads((ROOT / "MANIFEST.json").read_text())
    records = manifest["files"]
    listed = {record["path"] for record in records}
    require(len(listed) == len(records), "MANIFEST.json lists a path more than once")

    def is_local_worktree(rel: str) -> bool:
        name = Path(rel).name
        return (
            name in LOCAL_WORKTREE_NAMES
            or any(rel.startswith(prefix) for prefix in LOCAL_WORKTREE_PREFIXES)
        )

    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.name not in EXCLUDED_NAMES
        and not EXCLUDED_PARTS.intersection(path.parts) and path.suffix != ".pyc"
    }
    capsule_actual = {path for path in actual if not is_local_worktree(path)}
    require(listed == capsule_actual, (f"manifest coverage mismatch: "
                                       f"missing={sorted(capsule_actual - listed)}, "
                                       f"extra={sorted(listed - capsule_actual)}"))
    for record in records:
        path = ROOT / record["path"]
        require(path.is_file(), f"MANIFEST.json entry missing on disk: {record['path']}")
        require(record["bytes"] == path.stat().st_size,
                f"{record['path']} is {path.stat().st_size} bytes, manifest says {record['bytes']}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"],
                f"{record['path']} SHA-256 does not match MANIFEST.json")

    check_figure_contract()
    check_no_private_paths()
    version = manifest.get("version", "unknown")
    print(f"PASS: scReg-Eval audit capsule {version} artifacts, figures, manifest, and privacy are consistent")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValidationError, OSError, KeyError, ValueError) as error:
        print(f"FAIL: {type(error).__name__}: {error}", file=sys.stderr)
        sys.exit(1)
