# scReg-Eval red-team response workplan

**Status:** W0 done; W1: 0x10+0x13 done; remaining 0x11–0x12, 0x14–0x16  
**Inputs:** four independent red-teams (construct / stats / repro / PeerJ CS) + canvas `screg-eval-redteam`  
**Tip after merge fix:** `main` @ salvage-restored manuscript (`f75e70a`)  
**Goal:** strengthen the paper with text, figures, tables, and targeted experiments so reviewer attacks become design features / bounded claims.

Job IDs use hex codes (`0x01`…) and are mirrored in `.cursor/jobs.json` for agent settings.

---

## Priority waves

| Wave | Theme | Jobs |
|------|--------|------|
| W0 | Claim hygiene (text-only, no new compute) | `0x01`–`0x05` |
| W1 | Tables / figures from existing `results/*.public.json` | `0x10`–`0x16` |
| W2 | Small new analyses (reuse panels/null machinery) | `0x20`–`0x25` |
| W3 | Capsule / protocol integrity | `0x30`–`0x34` |
| W4 | Optional heavy / external (ENCODE, new panels) | `0x40`–`0x42` |

---

## Wave 0 — Claim hygiene (must land before resubmit)

### `0x01` Abstract reframe: Support ≠ recovery
- **Concerns:** C1/C3; construct #2/#5; stats #1/#2; PeerJ #2/#3
- **Work:** Rewrite Abstract Results so the first result sentence jointly states (i) full-spec dual-null count under disclosed degree hazard, (ii) non-degree dual-null = 0/13, (iii) positive dual-null rows fail usability and are not regulatory recovery.
- **Deliverable:** `paper/manuscript.tex` Abstract; sync PeerJ flat/source mirrors.
- **Acceptance:** No standalone “seven of thirteen … supported” without non-degree + usability in the same paragraph.

### `0x02` Kill “free of RNA covariation” overclaim
- **Concerns:** C2; construct #1/#3
- **Work:** Global replace framing → “edge weights built without expression statistics; panel and ATAC inputs remain expression-adjacent; proxy not claimed RNA-structure-independent.” Delete Intro “cannot be satisfied by co-expression.”
- **Deliverable:** Intro / Conclusion / contribution bullets.
- **Acceptance:** Zero remaining “free of RNA/expression covariation” absolute claims.

### `0x03` Estimand statement block
- **Concerns:** PeerJ #2; stats #5
- **Work:** New Methods subsection *Estimand*: fixed-panel unusualness under named confound × two named randomizations; not causal TF–target effect; not model ranking.
- **Deliverable:** ~½ page Methods; 2–3 sentence pointer in Discussion.
- **Acceptance:** Reader can quote a single paragraph answering “what question is answered?”

### `0x04` Rename usability → concordance / attention-likeness
- **Concerns:** construct #6; stats #8
- **Work:** Rename check in text + Fig. 4 title/caption; state it blocks analogy to attention-GRN stories, does **not** adjudicate residual regulatory content; fix “passers are null” vs brain Geneformer attention (+usability, negative dual-null).
- **Deliverable:** manuscript + `make_figs.R` labels for Fig. 4.
- **Acceptance:** Caption uses new name; contradictory passer called out in Results.

### `0x05` Bidirectional scope sentence + CS contribution front
- **Concerns:** construct #7; PeerJ #1/#5
- **Work:** Every “does not support global regulatory-capability claim” pairs with “equally cannot support field-level negative conclusion.” Lead Intro/Conclusion with protocol/software/report-discipline deliverables; biology as application scene.
- **Deliverable:** Intro contribution list; Conclusion; optional title softening (decision gate).
- **Acceptance:** Contribution (1) is protocol artifact, not “usable regulatory truth.”

---

## Wave 1 — Figures & tables from existing artifacts (high leverage)

### `0x10` Table 1b / Fig.3+: Protocol-pass matrix
- **Concerns:** C3; stats #2; PeerJ #3/#5
- **Data:** `results/fixed_panel_audit_v2.public.json` + usability Spearmans already in manuscript/JSON
- **Work:** New table (or Table 1 columns): dual-null full | dual-null non-degree | usability pass | multi-readout same-sign | vs-baseline | **protocol-pass**. Almost all rows fail protocol-pass — that is the point.
- **Figure option:** binary heatmap (rows × gates) next to Fig. 3.
- **Acceptance:** At least one manuscript sentence: “0/13 (or k/13) protocol-pass under predeclared gates.”

### `0x11` Fig.6 upgrade: full vs non-degree paired forest
- **Concerns:** C1; stats #1; repro #5
- **Data:** `results/spec_sensitivity_v2.public.json` / audit JSON both specs
- **Work:** Paired ρ points with sign-flip markers; annotate “dual-null both only under full.”
- **Acceptance:** Visual proves 7→0 dual-null without reading prose.

### `0x12` Effect-size & baseline panel (new Fig or Fig.3 inset)
- **Concerns:** C4; stats #3/#6; PeerJ #3
- **Data:** primary ρ, `brain_coexp_baseline_null_v2.public.json`, `pbmc_coexp_baseline_null_v2.public.json`
- **Work:** Bars/points for FM partial ρ vs co-expression baseline on same edge set; secondary axis or annotation for \(r^2=\rho^2\); note MC floor rate.
- **Acceptance:** Abstract can point to “same order as / below baseline” figure panel.

### `0x13` Finish table enrichment (Table 2 additive / Table 3 ranges)
- **Concerns:** PeerJ increment; prior `table-enrichment-plan.md`
- **Data:** `cross_tissue_additive_decomp_v2.public.json`, pertype summaries
- **Work:** Complete locked columns in enrichment plan; caption: cross-tissue ρ = motif/proximity stability, not regulatory truth (construct #4).
- **Acceptance:** Table 2 shows additive fraction 69–78% class numbers; Table 3 in manuscript.

### `0x14` Related-work comparison table (new Table)
- **Concerns:** PeerJ #4
- **Work:** 3–5 row table: scReg-Eval vs Kendiukhov scFM/SAE vs GeneRNIB vs BioLLM — reference has expression? degree control? null semantics? multi-readout? capsule?
- **Acceptance:** Intro cites table; “complement not compete” one-liner.

### `0x15` Scope card figure (new small schematic)
- **Concerns:** PeerJ #8; construct #4
- **Work:** One-panel “protocol instance scope”: ±2 kb × motif × ATAC; excludes distal/3D/ChIP; cell-type near-invariant.
- **Acceptance:** Limitations point to figure; reduces “reference too weak ⇒ models look empty” confusion.

### `0x16` Usability full scatter (Fig.4 fix)
- **Concerns:** stats #8
- **Data:** existing FM–coexp Spearmans
- **Work:** Plot all readouts including brain Geneformer attention +0.18 passer with dual-null negative; legend encodes dual-null status.
- **Acceptance:** No selective “passers are null” without that point visible.

---

## Wave 2 — Targeted new experiments / analyses

### `0x20` SAP freeze document + protocol-pass definition
- **Concerns:** stats #2; repro #8
- **Work:** Versioned `docs/SCREG_EVAL_SAP.md` (or strengthen `SCREG_EVAL_PROTOCOL.md` §5): Support ≠ capability; protocol-pass = dual-null ∧ usability ∧ non-degree same sign ∧ ρ > baseline (or explicit weaker rule). Hash-pin file in MANIFEST.
- **Experiment?** None — governance artifact.
- **Acceptance:** Manuscript cites SAP commit/hash; table columns match SAP.

### `0x21` FM−baseline superiority randomization (same confound)
- **Concerns:** stats #6
- **Work:** New small job: difference-of-partial-ρ (or paired edge residual contrast) under shared null draws; BH within tissue; report who beats baseline.
- **Deliverable:** `results/*_fm_vs_baseline_*.json` + table column / supplement.
- **Acceptance:** Multiplicity-controlled answer to “does any FM exceed co-expression baseline?”

### `0x22` Dual-null intersection operating characteristic
- **Concerns:** stats #4
- **Work:** Global-null simulation (or maxT/Westfall–Young) under shared proxy null batch; report empirical joint FDR of “both.”
- **Deliverable:** supplement note + optional Fig.5 inset.
- **Acceptance:** One calibrated number or “Support set unchanged under maxT.”

### `0x23` Estimand-matched null for non-degree
- **Concerns:** stats #5
- **Work:** For non-degree spec, add Mantel/label null comparison already emphasized; optionally a non-degree-preserving target shuffle; clarify dual rule only when nulls address same estimand.
- **Deliverable:** Methods caveat rewrite + optional sensitivity JSON.
- **Acceptance:** Discussion no longer hides only-D rejections behind “neither.”

### `0x24` Geometry / anisotropy diagnostic (lightweight)
- **Concerns:** construct #5 missing test
- **Work:** Per-TF embedding norm / cosine anisotropy vs degree; correlate with row-level contribution to partial ρ under full vs non-degree.
- **Deliverable:** supplement figure; one Results sentence if pattern holds.
- **Acceptance:** “anisotropy × degree” is no longer untested storytelling — or claim withdrawn.

### `0x25` Probe estimand cleanup + brain omission note
- **Concerns:** stats #7
- **Work:** Recompute q without random-init floor in family (sensitivity); state estimand Δρ vs ρ>MDE; explicit “brain not probed.”
- **Data:** `tf_probe_pair_*.public.json`
- **Acceptance:** Probe subsection cannot be skimmed as confirmatory FM>coexp.

---

## Wave 3 — Capsule / reproducibility (submission blockers)

### `0x30` Methods ↔ code: peak-count 40-cap
- **Concerns:** repro #1
- **Work:** Prefer **text fix** (peak-count = all peaks; GC uses ≤40) unless scientific reason to change code+rerun.
- **Acceptance:** `manuscript.tex` matches `src/fixed_panel_audit.py`.

### `0x31` Refresh MANIFEST + validate_artifacts PASS
- **Concerns:** repro #2
- **Work:** Rebuild PeerJ package from single SoT; update MANIFEST hashes; `python validate_artifacts.py` exits 0.
- **Acceptance:** README one-command gate green on clean tree.

### `0x32` FULL_RERUN boundary rewrite
- **Concerns:** repro #3/#4
- **Work:** Capsule `docs/FULL_RERUN.md` states heavy monorepo + external NPZ required; list SHA pins; no pretend `src/v2` in closed capsule unless vendored.
- **Acceptance:** Fresh clone instructions fail closed with clear missing-input message, not silent wrong path.

### `0x33` Panel ATAC-filter disclosure
- **Concerns:** repro #6; construct panel leakage
- **Work:** Methods sentence: universe requires ≥1 brain ATAC peak in window.
- **Acceptance:** “preregistered panel” not marketed as ATAC-orthogonal.

### `0x34` Soften or externally lock “preregistered”
- **Concerns:** repro #8
- **Work:** Either OSF/AsPredicted upload of panel+SAP, or wording → “protocol-frozen in-repo (hash …).”
- **Acceptance:** No overclaim of external preregistration without link.

---

## Wave 4 — Optional heavy (only if W0–W3 insufficient)

### `0x40` ENCODE / ChIP orthogonal calibration (if approvals clear)
- **Concerns:** construct missing calibration; PeerJ #8
- **Work:** Binding overlap enrichment for proxy edges; approval-gated.
- **Acceptance:** Pass/fail bound on capability language — not required for CS protocol paper if scope card (`0x15`) is strong.

### `0x41` No-ATAC-filter sensitivity panel
- **Concerns:** repro #6
- **Work:** Alternate gene universe without peak-presence filter; rerun primary audit subset.
- **Acceptance:** Dual-null qualitative stability statement.

### `0x42` allow_pickle purge on graph loaders
- **Concerns:** repro #7
- **Work:** Engineering sweep in monorepo `src/v2` loaders; capsule already mostly False.
- **Acceptance:** grep clean for `allow_pickle=True` on load paths used by release recipes.

---

## Suggested figure/table map (new or upgraded)

| ID | Asset | Wave | Primary concern |
|----|--------|------|-----------------|
| Fig.3+ / Table 1b | Protocol-pass heatmap/matrix | W1 | Support overclaim |
| Fig.6′ | Full vs non-degree paired forest | W1 | Degree hazard |
| Fig.3b / new | FM ρ vs baseline + r² | W1 | Tiny ρ |
| Fig.4′ | Full concordance scatter | W1 | Selective usability |
| Fig.S_scope | Scope card | W1 | Proxy narrowness |
| Table 2′ | Additive decomp columns | W1 | Construct stability |
| Table 3 | Per-type ranges | W1 | Enrichment plan |
| Table 4 | Related-work compare | W1 | PeerJ increment |
| Supp | FM>baseline test; dual-null OC | W2 | Stats #4/#6 |

---

## Execution order (safe default)

1. W0 text (`0x01`–`0x05`) on current salvage tip — unblocks honest visual QA.  
2. W1 tables/figs from public JSON (`0x10`–`0x16`) — no GPU.  
3. W3 capsule gate (`0x30`–`0x34`) before any PeerJ upload rebuild.  
4. W2 analyses (`0x20`–`0x25`) as supplement muscle.  
5. W4 only with explicit go-ahead (ENCODE / pickle / new panel).

## Non-goals

- Ranking models as better/worse regulators.  
- Claiming ENCODE-level ground truth without `0x40`.  
- Reintroducing retired bootstrap/MDE artifacts listed in `table-enrichment-plan.md` Forbidden.


## Execution log

- 2026-08-09: W0 (`0x01`–`0x05`) applied on `main` after salvage merge tip; local work branches deleted (fully merged).
- 2026-08-09: `0x10` protocol-pass matrix Table~\ref{tab:protocol} + Fig.~\ref{fig:protocol}; **0/13** protocol-pass under predeclared gates.
