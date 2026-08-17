#!/usr/bin/env python3
"""Build the sanitized public audit capsule from canonical project artifacts.

Copies an explicit allowlist from the canonical project tree, scrubs absolute
paths to placeholders, writes licenses and citation metadata, regenerates the
SHA-256 manifest, validates invariants, and emits a deterministic tarball.
"""
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERSION = "v0.4.0"
RELEASE_DATE = "2026-08-17"
CAPSULE_NAME = f"scfm-reg-audit-audit-capsule-{VERSION}"
RELEASE_DIR = ROOT / "release_candidate"
CAPSULE = RELEASE_DIR / CAPSULE_NAME
BRIDGE = ROOT / "release_private" / "ORIGINAL_TO_PUBLIC_HASH_BRIDGE.json"

# Order must match first-to-last \\input{figs/...} in paper/manuscript.tex.
CURRENT_FIGURES = (
    "fig_study_design.tex",
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
    "fig10_coverage_qc.tex",
)
CURRENT_TABLES = (
    "table5_related_work.tex",
    "table2_cross_tissue_observed.tex",
    "table1_primary_fixed_panel.tex",
    "table4_protocol_pass.tex",
    "table6_fm_vs_baseline.tex",
    "table3_pertype_ranges.tex",
    "table7_nondegree_null_pattern.tex",
    "table8_tf_probe_numeric.tex",
)
CURRENT_FRAGMENTS = CURRENT_FIGURES + CURRENT_TABLES
FIGURE_INPUT_RE = re.compile(r"\\input\{figs/([^}]+\.tex)\}")
WRAPPER_INPUT_RE = re.compile(r"\\input\{([^}/]+\.tex)\}")
FLAT_FIGURE_RE = re.compile(r"\\includegraphics\[[^]]*\]\{(Figure\d+\.pdf)\}")
FLAT_APPENDIX_RE = re.compile(r"\\includegraphics\[[^]]*\]\{(FigureA\d+\.pdf)\}")

# (source relative to ROOT, destination relative to capsule)
COPY_FILES = [
    ("README.md", "README_PROJECT.md"),
    ("LICENSE", "LICENSE"),
    ("LICENSE-CONTENT.md", "LICENSE-CONTENT.md"),
    ("LICENSING.md", "LICENSING.md"),
    ("data/manifest/shared_genes.v2.json", "data/manifest/shared_genes.v2.json"),
    ("docs/PAPER_OUTLINE.md", "docs/PAPER_OUTLINE.md"),
    ("results/v2/LEGACY_INFERENCE_NOTE.md", "docs/LEGACY_INFERENCE_NOTE.md"),
    ("docs/SCREG_EVAL_PROTOCOL.md", "docs/SCREG_EVAL_PROTOCOL.md"),
    ("docs/TF_PROBE_RESULT.md", "docs/TF_PROBE_RESULT.md"),
    ("docs/FULL_RERUN.md", "docs/FULL_RERUN.md"),
    ("docs/FULL_RUN_LOG.md", "docs/FULL_RUN_LOG.md"),
    ("docs/NOTICE.md", "docs/NOTICE.md"),
    ("paper/NOTICE.md", "paper/NOTICE.md"),
    ("results/v2/NOTICE.md", "results/NOTICE.md"),
    ("ENVIRONMENT.example", "ENVIRONMENT.example"),
    ("paper/manuscript.tex", "paper/manuscript.tex"),
    ("paper/references.bib", "paper/references.bib"),
    ("paper/make_figs.R", "paper/make_figs.R"),
    ("paper/make_panel_data.py", "paper/make_panel_data.py"),
    ("paper/panel_data.json", "paper/panel_data.json"),
    ("paper/figs_preview.tex", "paper/figs_preview.tex"),
    ("REPO_TOPOLOGY.md", "REPO_TOPOLOGY.md"),
    ("docs/SCREG_EVAL_SAP.md", "docs/SCREG_EVAL_SAP.md"),
    ("paper/CANONICAL_BUILD.md", "paper/CANONICAL_BUILD.md"),
    ("src/emit_wave5_tables.py", "src/emit_wave5_tables.py"),
    ("src/w2_public_analyses.py", "src/w2_public_analyses.py"),
    ("src/v2/tests/test_coexp_baseline_null.py", "src/tests/test_coexp_baseline_null.py"),
    ("src/v2/tests/test_pair_probe.py", "src/tests/test_pair_probe.py"),
    ("src/v2/tests/test_pbmc_cache.py", "src/tests/test_pbmc_cache.py"),
    ("src/tests/test_validate_artifacts.py", "src/tests/test_validate_artifacts.py"),
    ("src/v2/fixed_panel_audit.py", "src/fixed_panel_audit.py"),
    ("src/v2/run_fixed_panel_audit.py", "src/run_fixed_panel_audit.py"),
    ("src/v2/pbmc_cache.py", "src/pbmc_cache.py"),
    ("src/v2/benchmark_n99.py", "src/benchmark_n99.py"),
    ("src/v2/pbmc_uce_eval_v2.py", "src/pbmc_uce_eval_v2.py"),
    ("src/v2/brain_coexp_baseline_null.py", "src/brain_coexp_baseline_null.py"),
    ("src/v2/pbmc_coexp_baseline_null.py", "src/pbmc_coexp_baseline_null.py"),
    ("src/v2/pair_probe_stats.py", "src/pair_probe_stats.py"),
    ("src/v2/run_pair_probe.py", "src/run_pair_probe.py"),
    ("src/v2/tests/test_fixed_panel_audit.py", "src/tests/test_fixed_panel_audit.py"),
]
COPY_FILES += [
    (f"paper/figs/{name}", f"paper/figs/{name}")
    for name in CURRENT_FRAGMENTS
]

PUBLIC_JSONS = [
    "fixed_panel_audit_v2.json",
    "fixed_panel_signal_injection_v2.json",
    "inference_status_v2.json",
    "brain_coexp_baseline_null_v2.json",
    "pbmc_coexp_baseline_null_v2.json",
    "tf_probe_pair_stats_v2.json",
    "tf_probe_pair_sensitivity_v2.json",
    "proxy_celltype_invariance_v2.json",
    "injection_subdivided_v2.json",
    "cross_tissue_additive_decomp_v2.json",
    # figure-atlas inputs (v0.2.6+): every JSON make_figs.R reads
    "tf_probe_pair_eval_v2.json",
    "spec_sensitivity_v2.json",
    "marginal_vs_adjusted_v2.json",
    "effect_vs_injection_scale_v2.json",
    "degree_preserving_null_v2.json",
    "readout_attention_v2.json",
    "insilico_ko_v2.json",
    "verify_brain_attention_omission_v2.json",
    "pertype_fm_v2.json",
]

# Already-public derivatives retained from the v0.3.0 tree (no private raw sibling).
EXTRA_PUBLIC_JSONS = [
    "dual_null_oc_independence_v2.public.json",
    "encode_proxy_calibration_v1.public.json",
    "fm_vs_baseline_observed_v2.public.json",
    "fm_vs_baseline_shared_null_v2.public.json",
    "no_atac_filter_manifest_sensitivity_v1.public.json",
    "nondegree_null_pattern_v2.public.json",
    "tf_probe_contrasts_no_floor_v2.public.json",
]

PATH_REPLACEMENTS = [
    (str(ROOT), "${SCFM_PROJECT_ROOT}"),
    ("/home/zeyufu/Desktop/data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad",
     "${SCFM_BRAIN_ATAC}"),
    ("/home/zeyufu/Desktop/data", "${SCFM_DATA_ROOT}"),
    ("/home/zeyufu", "${HOME}"),
]

TEXT_SUFFIXES = {".json", ".md", ".tex", ".py", ".R", ".bib", ".cff", ".txt"}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sanitize_text(text):
    for old, new in PATH_REPLACEMENTS:
        text = text.replace(old, new)
    return text


def copy_with_sanitization(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix in TEXT_SUFFIXES:
        dst.write_text(sanitize_text(src.read_text()))
    else:
        shutil.copy2(src, dst)


def figure_inputs(path):
    return tuple(FIGURE_INPUT_RE.findall(path.read_text()))


def validate_figure_contract(root):
    manuscript_inputs = figure_inputs(root / "paper/manuscript.tex")
    manuscript_figures = tuple(name for name in manuscript_inputs if name.startswith("fig"))
    manuscript_tables = tuple(name for name in manuscript_inputs if name.startswith("table"))
    preview_inputs = figure_inputs(root / "paper/figs_preview.tex")

    assert manuscript_figures == CURRENT_FIGURES, (
        f"manuscript figure order mismatch: {manuscript_figures}")
    assert manuscript_tables == CURRENT_TABLES, (
        f"manuscript table order mismatch: {manuscript_tables}")
    assert preview_inputs == CURRENT_FRAGMENTS, (
        f"preview/manuscript fragment mismatch: {preview_inputs}")

    wrappers = root / "paper/submission_peerj/internal/figure_build"
    wrapper_inputs = tuple(
        WRAPPER_INPUT_RE.findall((wrappers / f"Figure{index}.tex").read_text())[0]
        for index in range(1, len(CURRENT_FIGURES) + 1)
    )
    assert wrapper_inputs == CURRENT_FIGURES, (
        f"PeerJ wrapper figure order mismatch: {wrapper_inputs}")

    flat_tex = root / "paper/submission_peerj/flat_upload/manuscript.tex"
    flat_figures = tuple(FLAT_FIGURE_RE.findall(flat_tex.read_text()))
    expected_flat = tuple(f"Figure{index}.pdf" for index in range(1, len(CURRENT_FIGURES) + 1))
    assert flat_figures == expected_flat, (
        f"PeerJ flat figure order mismatch: {flat_figures}")
    flat_appendix = tuple(FLAT_APPENDIX_RE.findall(flat_tex.read_text()))
    assert flat_appendix == ("FigureA1.pdf", "FigureA2.pdf", "FigureA3.pdf"), (
        f"PeerJ flat appendix figure mismatch: {flat_appendix}")


def write_citation():
    citation = f"""cff-version: 1.2.0
title: "scReg-Eval: a fixed-panel audit of regulatory alignment in single-cell RNA foundation-model gene graphs"
message: "If you use this audit capsule, please cite it as below."
type: software
authors:
  - family-names: Fu
    given-names: Zeyu
    orcid: "https://orcid.org/0009-0001-8329-0108"
    affiliation: "Army Medical University"
    email: "fuzeyu99@126.com"
license: "MIT AND CC-BY-4.0"
notes: "Top-level MIT applies to original software code only. Manuscript, figures, tables, and derived public results are CC BY 4.0; see LICENSING.md."
keywords:
  - single-cell RNA-seq
  - foundation models
  - gene regulatory networks
  - co-expression confounding
  - randomization inference
doi: "10.5281/zenodo.21724336"
url: "https://github.com/PeterPonyu/scfm-reg-audit"
version: "{VERSION.lstrip('v')}"
date-released: "{RELEASE_DATE}"
"""
    (CAPSULE / "CITATION.cff").write_text(citation)


def write_readme():
    readme = f"""# scReg-Eval fixed-panel audit capsule {VERSION}

This is a sanitized **audit capsule** for the scReg-Eval manuscript. It validates the published
numerical artifacts and carries the audited statistical code; it is not a raw-data reproduction
environment.

## Included

- fixed-panel statistical implementation, PBMC UCE pipeline, baseline and probe scripts, and tests;
- frozen 446-TF / 1,200-gene manifest;
- path-scrubbed authoritative JSON derivatives (pooled audit, injection, inference status,
  co-expression baselines, TF-disjoint probe, cell-type invariance, subdivided calibration);
- manuscript source, active data-driven figure fragments, and the figure generator;
- licenses (MIT for original code, CC BY 4.0 for manuscript/figures/derived results), citation
  metadata, and this validator with a SHA-256 manifest.

## Excluded

Raw/processed H5AD/H5 data, genome FASTA, model weights, cached NPZ graphs, model-inference/vendor
code, pilot code, session state, build logs/caches, internal review records, and the four retired
numerical JSON files are excluded. See `docs/PAPER_OUTLINE.md` and `docs/LEGACY_INFERENCE_NOTE.md`.

## Validate

```bash
python validate_artifacts.py
```

## Test

```bash
pip install numpy scipy anndata pyfaidx  # test dependencies
python -m unittest src.tests.test_fixed_panel_audit
```

The bundled suite is the statistical-contract suite (`test_fixed_panel_audit.py`). Legacy-hash,
model-scope, and real-cache integration cases skip automatically when their retired or
machine-local inputs are absent; the development repository runs the full suite including
model/cache contract tests. The release builder runs this command after every rebuild and fails
on any error.

## Reproduction boundary

The capsule validates published numerical artifacts. Re-running production requires external public
datasets and cached model/readout graphs not redistributed here. See `docs/FULL_RERUN.md` in the
development repository for the full-rerun recipe. Set `SCFM_BRAIN_ATAC` to the brain
ATAC H5AD path. Provenance strings use `${{SCFM_PROJECT_ROOT}}`, `${{SCFM_DATA_ROOT}}`, and
`${{SCFM_BRAIN_ATAC}}` placeholders.
"""
    (CAPSULE / "README.md").write_text(readme)


def install_validator():
    # Tip SoT is repo-root validate_artifacts.py; keep src/v2/validate_capsule.py
    # as a synced mirror for older build scripts / offline copies.
    tip = ROOT / "validate_artifacts.py"
    mirror = Path(__file__).resolve().parent / "validate_capsule.py"
    src = tip if tip.is_file() else mirror
    shutil.copy2(src, CAPSULE / "validate_artifacts.py")
    if tip.is_file() and mirror.is_file():
        tip_text = tip.read_text()
        if mirror.read_text() != tip_text:
            mirror.write_text(tip_text)


def write_manifest():
    records = []
    for path in sorted(CAPSULE.rglob("*")):
        if (path.is_file() and path.name != "MANIFEST.json"
                and "__pycache__" not in path.parts and path.suffix != ".pyc"):
            records.append({
                "path": path.relative_to(CAPSULE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    # MANIFEST.json is content-addressed only: no build-time metadata inside the
    # tarball, so the same content yields the same tarball hash at any commit.
    (CAPSULE / "MANIFEST.json").write_text(
        json.dumps({"capsule": CAPSULE_NAME, "version": VERSION, "files": records}, indent=1) + "\n")
    builder_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import platform
    (RELEASE_DIR / f"{CAPSULE_NAME}.BUILD_INFO.json").write_text(json.dumps({
        "builder": {
            "script": "src/v2/build_release_capsule.py",
            "commit": builder_commit,
            "python": platform.python_version(),
        },
        "tarball": f"{CAPSULE_NAME}.tar.gz",
    }, indent=1) + "\n")


def write_bridge():
    records = []
    for name in PUBLIC_JSONS:
        original = ROOT / "results/v2" / name
        public = CAPSULE / "results" / (name.replace(".json", ".public.json"))
        records.append({
            "original": f"results/v2/{name}",
            "original_sha256": sha256_file(original),
            "public": f"results/{public.name}",
            "public_sha256": sha256_file(public),
            "transform": "absolute path scrub to placeholders",
            "release": VERSION,
        })
    BRIDGE.parent.mkdir(parents=True, exist_ok=True)
    # Canonical JSON: sorted keys, compact separators, trailing newline. The validator
    # and future tooling must reproduce this byte representation exactly.
    canonical = json.dumps({"release": VERSION, "records": records},
                           sort_keys=True, separators=(",", ":")) + "\n"
    BRIDGE.write_text(canonical)


def deterministic_tarball():
    archive = RELEASE_DIR / f"{CAPSULE_NAME}.tar.gz"
    with open(archive, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar:
                for path in sorted(CAPSULE.rglob("*")):
                    if not path.is_file():
                        continue
                    info = tar.gettarinfo(
                        str(path),
                        arcname=f"{CAPSULE_NAME}/{path.relative_to(CAPSULE).as_posix()}")
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = 0
                    info.mode = 0o644
                    with open(path, "rb") as handle:
                        tar.addfile(info, handle)
    digest = sha256_file(archive)
    (RELEASE_DIR / f"{CAPSULE_NAME}.tar.gz.sha256").write_text(f"{digest}  {archive.name}\n")
    return archive, digest


def run_capsule_validator():
    result = subprocess.run(
        ["python3", "validate_artifacts.py"],
        cwd=CAPSULE, capture_output=True, text=True)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError("capsule validator failed:\n" + detail)
    print("capsule validator:", result.stdout.strip())


def run_capsule_tests():
    result = subprocess.run(
        ["python3", "-m", "unittest", "src.tests.test_fixed_panel_audit"],
        cwd=CAPSULE, capture_output=True, text=True)
    for cache_dir in CAPSULE.rglob("__pycache__"):
        shutil.rmtree(cache_dir)
    tail = (result.stderr or result.stdout).strip().splitlines()[-4:]
    if result.returncode != 0:
        raise RuntimeError("capsule test suite failed:\n" + "\n".join(tail))
    print("capsule tests:", " | ".join(tail[-2:]))


def main():
    validate_figure_contract(ROOT)
    if CAPSULE.exists():
        shutil.rmtree(CAPSULE)
    CAPSULE.mkdir(parents=True)

    for src_rel, dst_rel in COPY_FILES:
        src = ROOT / src_rel
        if not src.exists():
            raise FileNotFoundError(src_rel)
        copy_with_sanitization(src, CAPSULE / dst_rel)

    for name in PUBLIC_JSONS:
        src = ROOT / "results/v2" / name
        if not src.exists():
            raise FileNotFoundError(name)
        dst = CAPSULE / "results" / name.replace(".json", ".public.json")
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(sanitize_text(src.read_text()))

    for name in EXTRA_PUBLIC_JSONS:
        src = ROOT / "results" / name
        if not src.exists():
            raise FileNotFoundError(name)
        dst = CAPSULE / "results" / name
        dst.write_text(sanitize_text(src.read_text()))

    write_citation()
    write_readme()
    install_validator()
    run_capsule_tests()
    write_manifest()
    run_capsule_validator()
    write_bridge()
    archive, digest = deterministic_tarball()
    print(f"capsule: {CAPSULE}")
    print(f"archive: {archive} sha256={digest}")


if __name__ == "__main__":
    main()
