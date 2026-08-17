# scReg-Eval — Statistical Analysis Plan (SAP)

**Status:** protocol-frozen in-repo (Wave 2 / job `0x20`)  
**Version:** SAP-v1.0 (2026-08-09)  
**Panel pins (unchanged from audit v2):**

| pin | SHA-256 |
|-----|---------|
| `manifest_sha256` | `6b203fcfab45dc600f84d2149c7f5f94e1a876f584529a0b465694e170b4f848` |
| `tf_panel_sha256` | `b07ae73888cd2e075cd1992f73b5fac9a2fa9c5d8f73c596eac653d259eae8da` |

This SAP freezes *definitions and reporting rules* for the fixed-panel audit. It is
**not** an external OSF/AsPredicted preregistration. It complements
`docs/SCREG_EVAL_PROTOCOL.md` (construct) and `docs/FULL_RERUN.md` (compute boundary).

---

## 1. Estimand

**Fixed-panel unusualness:** on the protocol-frozen 446-TF × 1,200-gene panel, is the
partial-Spearman alignment between a frozen foundation-model gene graph and the
accessibility/motif proxy unusually large under (i) a named confound design and
(ii) two named Monte Carlo randomizations, after within-family BH control?

Dual-null **Support** answers that panel-local question only. It is **not**:

- causal TF–target recovery,
- a ranking of foundation models,
- residual regulatory content after every possible confound,
- a field-level negative about scFMs.

---

## 2. Primary statistic and confound specs

- **Statistic:** partial Spearman via residualized ranks (Pearson of residuals).
- **Full design:** intercept + ranked co-expression + z(peakcount) + z(genelen) +
  z(detection) + z(GC) + z(TF out-degree) + z(target in-degree).
  - Peak count = **all** linked peaks in the gene window.
  - GC = mean over **at most 40** peaks/gene (sampling only).
  - Degree columns are **proxy-derived** (outcome functions) → validity hazard.
- **Non-degree design (prespecified sensitivity):** full minus the two degree columns.
- Matching column sets for observed and randomized statistics.

---

## 3. Nulls and multiplicity

| Null | What is randomized | Plus-one \(N=999\) |
|------|--------------------|--------------------|
| Gene-label (Mantel-style) | Proxy edge endpoints | \(p_{\mathrm{MC}}=(1+\#\{|T_b|\ge|T_{\mathrm{obs}}|\})/1000\) |
| Degree-preserving row shuffle | Within-TF target assignment | same |

- One shared proxy randomization batch per replicate across model rows in a family.
- **BH families (8):** tissue × confound-spec × null.
- **Dual-null Support:** both \(q_M<0.05\) and \(q_D<0.05\) in the same confound-spec.
- Opposite-sign dual support reported as `both (neg)`.
- Co-expression baseline is **outside** FM BH families (tested on the same edge set
  with the co-expression column removed from the design).

### Estimand-matched dual rule (non-degree)

Under non-degree, the row-shuffle null still preserves TF out-degree while the
statistic does **not** condition on degree. **D-only** rejections under non-degree
are therefore not dual-Support and must not be sold as “Support under the non-degree
sensitivity.” Dual-Support requires **both** nulls to reject under the same named design.

---

## 4. Protocol-pass (reporting gate, Table 4 / Fig. 12)

Protocol-pass is a **predeclared conjunction** for each of the 13 FM/readout rows.
It is stricter than dual-null Support and is the gate used when claiming
“recovery-like” alignment.

| Gate | Definition |
|------|------------|
| Dual full | \(q_M<0.05\) and \(q_D<0.05\) under **full** |
| Dual non-deg | same under **non-degree** (reported; typically all fail) |
| Concordance | FM–co-expression Spearman \(>0\) (attention-likeness) |
| ND same sign | \(\mathrm{sign}(\rho_{\mathrm{full}})=\mathrm{sign}(\rho_{\mathrm{nondeg}})\) |
| Multi-RO sign | full-spec \(\rho\) signs agree within model family × tissue |
| \(\rho>\) baseline | \(\rho_{\mathrm{full}} > \rho_{\mathrm{coexp,baseline}}\) for that tissue |
| **Protocol-pass** | Dual full ∧ Concordance ∧ ND same sign ∧ Multi-RO sign ∧ \(\rho>\) baseline |

**Interpretation:** dual-null Support **without** protocol-pass is unusual residual
geometry under disclosed degree hazard — **not** regulatory recovery.

Empirical freeze (audit v2 public JSON): **0/13** protocol-pass.

---

## 5. FM vs co-expression baseline (job `0x21`)

- Observed \(\Delta\rho = \rho_{\mathrm{FM,full}} - \rho_{\mathrm{baseline,tissue}}\).
- Report per-row \(\Delta\rho\), `beats_baseline` (\(\Delta\rho>0\)), and intersection with
  dual-null Support.
- Shared-null randomization of \(\Delta\rho\) requires monorepo NPZ caches; the capsule
  reports **observed** superiority only, clearly labeled
  `method=observed_delta_v1` in `results/fm_vs_baseline_observed_v2.public.json`.

---

## 6. Dual-null operating characteristic (job `0x22`)

Under a **global null** simulation with independent Uniform\(p\) draws in each BH
family (size = number of FM rows in that tissue×full family), report the Monte Carlo
rate of ≥1 dual-null Support row. Shared-batch dependence is **not** simulated here;
independence is a reference OC, not a claim of exact FDR control under the audit’s
shared null stream. See `results/dual_null_oc_independence_v2.public.json`.
The manuscript does **not** print this independence OC as the Type I of Support;
Dual-null Support is a reporting gate, not an independence Type I rate.

---

## 7. Probe arm (job `0x25`)

- Estimand: supervised TF-disjoint recoverability of proxy profiles from edge features
  on **paired PBMC multiome only** — not the eight brain rows.
- Primary contrasts vs co-expression use per-TF sign-flip \(p\) with BH within the
  contrast family.
- Sensitivity: BH recomputed **excluding** `random_floor` from the contrast family
  (`results/tf_probe_contrasts_no_floor_v2.public.json`).
- Probe results are **not** confirmatory evidence that all FMs beat co-expression.
- Figure 10 gene-label \(q_M\) (when drawn) is **display-only**, not a BH family.
  The probe contrast family remains paired sign-flip \(q_{\mathrm{flip}}\).

---

## 8. Geometry / anisotropy (job `0x24`)

Per-TF embedding-norm / cosine anisotropy vs degree is **not** a frozen primary
analysis in this SAP. Manuscript language must not treat it as a tested mechanism
unless a protocol-frozen diagnostic JSON is added later.

---

## 9. Manuscript mapping

| Asset | SAP section |
|-------|-------------|
| Table 1 Support column | §3 Dual-null Support |
| Table 4 / Fig. 12 protocol-pass | §4 |
| Fig. 3D / baseline text | §5 |
| Non-degree Methods caveat | §3 estimand-matched dual rule |
| Probe subsection | §7 |

File hash of this SAP is pinned in `MANIFEST.json` after each freeze edit.
