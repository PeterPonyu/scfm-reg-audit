#!/usr/bin/env python3
"""Emit Wave-5 numeric table fragments from capsule-safe public JSON.

Writes:
  paper/figs/table6_fm_vs_baseline.tex   (observed Δ + shared-null q when available)
  paper/figs/table7_nondegree_null_pattern.tex
  paper/figs/table8_tf_probe_numeric.tex

Sources (tip results/):
  fm_vs_baseline_observed_v2.public.json
  fm_vs_baseline_shared_null_v2.public.json  (optional until Option B lands)
  nondegree_null_pattern_v2.public.json
  tf_probe_pair_stats_v2.public.json
  tf_probe_contrasts_no_floor_v2.public.json
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGS = ROOT / "paper" / "figs"

MODEL_LABEL = {
    "geneformer_embed": "Geneformer embed",
    "geneformer_attn": "Geneformer attention",
    "geneformer_ko_raw": "Geneformer KO",
    "geneformer_ko_posctrl": "Artifact-corrected KO",
    "scFoundation_encoder": "scFoundation encoder",
    "UCE_encoder": "UCE encoder",
    "scGPT_encoder": "scGPT encoder",
    "random_init_floor": "Random-init floor",
    "random_floor": "Random-init floor",
    "co_expression": "Co-expression",
}

PROBE_FAMILY_ORDER = (
    "co_expression",
    "geneformer_embed",
    "geneformer_attn",
    "scGPT_encoder",
    "UCE_encoder",
    "random_floor",
)

PATTERN_LABEL = {
    "dual": "dual",
    "M_only": "M-only",
    "D_only": "D-only",
    "neither": "neither",
}

TISSUE_LABEL = {
    "brain": "Brain",
    "pbmc": "PBMC",
}


def load_json(name: str) -> dict:
    return json.loads((RESULTS / name).read_text())


def fmt_signed(x: float, digits: int = 5) -> str:
    return f"${x:+.{digits}f}$"


def fmt_q(x: float | None, digits: int = 3) -> str:
    if x is None:
        return "--"
    return f"${x:.{digits}f}$"


def fmt_p(p: float | None) -> str:
    if p is None:
        return "--"
    if p <= 0.001:
        return "$<0.001$"
    return f"${p:.3f}$"


def yn(flag: bool) -> str:
    return "yes" if flag else "no"


def emit_table6(path: Path) -> None:
    """Emit FM−baseline table; prefer shared-null JSON for ρ/Δ when present."""
    obs = load_json("fm_vs_baseline_observed_v2.public.json")
    shared_path = RESULTS / "fm_vs_baseline_shared_null_v2.public.json"
    shared = json.loads(shared_path.read_text()) if shared_path.is_file() else None
    obs_by = {(r["tissue"], r["model_label"]): r for r in obs["rows"]}
    # Prefer shared-null row order/values so printed Δ matches the tested Δ.
    rows_src = shared["rows"] if shared is not None else obs["rows"]

    body = []
    for r in rows_src:
        key = (r["tissue"], r["model_label"])
        o = obs_by.get(key, {})
        if shared is not None:
            rho = float(r["rho_fm_full"])
            delta = float(r["delta_rho"])
            beats = bool(r["beats_baseline_observed"])
            dual = bool(o.get("dual_null_full", False))
            dual_beats = bool(dual and beats)
            q_cell = fmt_q(float(r["bh_q"]), 3)
            sig = yn(bool(r["shared_null_significant"]))
        else:
            rho = float(r["rho_full"])
            delta = float(r["delta_rho"])
            beats = bool(r["beats_baseline"])
            dual = bool(r["dual_null_full"])
            dual_beats = bool(r["dual_and_beats_baseline"])
            q_cell = "--"
            sig = "--"
        body.append(
            " & ".join(
                [
                    TISSUE_LABEL.get(r["tissue"], r["tissue"].title()),
                    MODEL_LABEL.get(r["model_label"], r["model_label"]),
                    fmt_signed(rho, 5),
                    fmt_signed(delta, 5),
                    yn(dual),
                    yn(beats),
                    yn(dual_beats),
                    q_cell,
                    sig,
                ]
            )
            + " \\\\"
        )

    lines = [
        "\\begin{tabular}{llrrccccc}",
        "\\toprule",
        (
            "Tissue & Readout & $\\rho_{\\mathrm{FM}}$ & $\\Delta\\rho$ vs base & "
            "Dual-null & Beats base & Dual$\\wedge$beats & "
            "$q_{\\Delta}$ (shared) & Shared sig \\\\"
        ),
        "\\midrule",
        *body,
        "\\bottomrule",
        "\\end{tabular}",
        "",
    ]
    path.write_text("\n".join(lines))


def emit_table7(path: Path) -> None:
    data = load_json("nondegree_null_pattern_v2.public.json")
    rows = []
    for r in data["rows"]:
        tissue = TISSUE_LABEL.get(r["tissue"], str(r["tissue"]).title())
        model = MODEL_LABEL.get(r["model_label"], r["model_label"])
        pattern = PATTERN_LABEL.get(r["pattern"], r["pattern"])
        rows.append(
            f"{tissue} & {model} & {fmt_signed(float(r['rho']))} & "
            f"{fmt_q(float(r['q_M']))} & {fmt_q(float(r['q_D']))} & {pattern} \\\\"
        )
    lines = [
        "\\begin{tabular}{llrrrl}",
        "\\toprule",
        "Tissue & Model & $\\rho$ & $q_M$ & $q_D$ & Pattern \\\\",
        "\\midrule",
        *rows,
        "\\bottomrule",
        "\\end{tabular}",
        "",
    ]
    path.write_text("\n".join(lines))


def emit_table8(path: Path) -> None:
    stats = load_json("tf_probe_pair_stats_v2.public.json")
    no_floor = load_json("tf_probe_contrasts_no_floor_v2.public.json")
    families = stats["families"]
    contrasts = stats["contrasts_vs_baseline"]
    nf_contrasts = no_floor["contrasts_without_random_floor"]

    rows = []
    for fam in PROBE_FAMILY_ORDER:
        if fam not in families:
            raise KeyError(f"missing probe family {fam!r} in pair stats")
        fam_stats = families[fam]
        label = MODEL_LABEL.get(fam, fam)
        rho = float(fam_stats["adjusted_rho_mean"])
        if fam == "co_expression":
            rows.append(
                f"{label} & {fmt_signed(rho)} & -- & "
                f"{fmt_p(float(fam_stats['mantel_p']))} & "
                f"{fmt_q(float(fam_stats['mantel_q']), digits=4)} & -- \\\\"
            )
            continue

        contr = contrasts[fam]
        delta = float(contr["paired_delta_mean"])
        p_flip = float(contr["signflip_p"])
        q_flip = float(contr["signflip_q"])
        if fam in nf_contrasts:
            q_nf = float(nf_contrasts[fam]["signflip_q_family_without_random_floor"])
            q_nf_tex = fmt_q(q_nf, digits=4)
        else:
            q_nf_tex = "--"
        # Prefer four decimals when q carries the US2 UCE 0.0325→0.052 story.
        q_flip_tex = fmt_q(q_flip, digits=4)
        rows.append(
            f"{label} & {fmt_signed(rho)} & {fmt_signed(delta)} & "
            f"{fmt_p(p_flip)} & {q_flip_tex} & {q_nf_tex} \\\\"
        )

    lines = [
        "\\begin{tabular}{lrrrrr}",
        "\\toprule",
        (
            "Family & $\\rho$ & $\\Delta\\rho$ vs base & "
            "$p_{\\mathrm{flip}}$ & $q_{\\mathrm{flip}}$ & "
            "$q$ (no floor) \\\\"
        ),
        "\\midrule",
        *rows,
        "\\bottomrule",
        "\\end{tabular}",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    out6 = FIGS / "table6_fm_vs_baseline.tex"
    out7 = FIGS / "table7_nondegree_null_pattern.tex"
    out8 = FIGS / "table8_tf_probe_numeric.tex"
    emit_table6(out6)
    emit_table7(out7)
    emit_table8(out8)
    print(f"wrote {out6.relative_to(ROOT)}")
    print(f"wrote {out7.relative_to(ROOT)}")
    print(f"wrote {out8.relative_to(ROOT)}")
    shared = RESULTS / "fm_vs_baseline_shared_null_v2.public.json"
    if shared.is_file():
        n = json.loads(shared.read_text()).get("n_perm")
        print(f"table6 includes shared-null columns (n_perm={n})")
    else:
        print("table6 shared-null columns are placeholders (--); JSON missing")


if __name__ == "__main__":
    main()
