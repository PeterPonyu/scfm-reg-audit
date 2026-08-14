# Design

Durable design contract for the **GitHub Pages materials site** of this academic audit
capsule. Implementation must cite this file. This document does **not** implement the
site, and it does **not** change manuscript or figure-generator sources.

A product-requirements note may land separately at `.omx/plans/prd-github-pages.md`.
If that file exists, do not overwrite it. This `DESIGN.md` owns brand, IA, visual
language, components, accessibility, and content voice.

## Source of truth

- Status: Draft
- Last refreshed: 2026-08-13
- Primary product surfaces:
  - Planned GitHub Pages site for paper **materials / products / resources**
    (expected project URL `https://peterponyu.github.io/scfm-reg-audit/`).
  - Two independent reading paths for the same scReg-Eval audit: **PeerJ Computer
    Science** and **Frontiers in Genetics**.
  - Shared figure stems, review PDFs, capsule validation, and code/data pointers.
  - Not a second manuscript. Not a dashboard. Not a marketing landing page.
- Evidence reviewed:
  - **Observed (no existing web UI):** no repo-root `DESIGN.md` before this draft;
    no `docs/design*`, `docs/ux*`, `docs/frontend*`; no `package.json`, `index.html`,
    `_config.yml`, `.github/` workflows, `CNAME`, or `gh-pages` / `github.io` traces.
    PNG files under `paper/figs_extension/visual_qa/` and
    `docs/reports/extension-plans/visual-qa-human/` are local TikZ visual-QA
    scratch, not a public gallery.
  - **Capsule / product:** `README.md` (v0.3.0 audit capsule; validates artifacts;
    not a raw-data reproduction environment); `CITATION.cff` (title, Zenodo
    `10.5281/zenodo.21724336`, GitHub `PeterPonyu/scfm-reg-audit`, MIT AND CC-BY-4.0);
    `LICENSING.md`, `LICENSE`, `LICENSE-CONTENT.md`, `docs/NOTICE.md`.
  - **Two report paths:** `paper/CANONICAL_BUILD.md`;
    `paper/submission_frontiers_genetics/README.md`; `paper/submission_peerj/README.md`;
    `paper/submission_peerj/HUMAN_GATES.md`. PeerJ SoT is `paper/manuscript.tex` →
    `paper/manuscript.pdf`. Frontiers SoT is
    `paper/submission_frontiers_genetics/manuscript.tex` → `manuscript.pdf` in that
    folder. Frontiers is a **full conversion** (main Figures 1–13 + Tables 1–8 +
    Appendix A1–A3), not a three-figure-only paper.
  - **Claims / protocol:** `docs/PAPER_OUTLINE.md` (forbidden claim patterns);
    `docs/SCREG_EVAL_PROTOCOL.md`; `docs/SCREG_EVAL_SAP.md` (protocol-pass
    conjunction; dual-null Support ≠ recovery); `docs/FULL_RERUN.md` (capsule
    validate-only vs experiment-workspace full rerun); `docs/LEGACY_INFERENCE_NOTE.md`
    (retired bootstrap/MDE — do not surface as current results).
  - **Figure language:** `paper/make_figs.R` tokens `AQUA <- "#1b8f75"`,
    `BLUE <- "#2a78d6"`, `YELLOW <- "#d99a00"`, `VIOLET <- "#5946b2"`,
    `RED <- "#d84a4a"`, `MUTED <- "grey45"`, Fig 12 fail fill `grey88`, conjunction
    border `#5a6570`, header ink `#1a1a1a`. `paper/figs/fig12_protocol_pass_matrix.tex`
    encodes pass `RGB{27,143,117}` / fail `gray{0.88}`. Frontiers caption: “Aqua pass /
    grey fail”. Protocol-pass frozen **0/13**.
  - **Typography in PDFs (not the web stack):** PeerJ / TikZ contract uses
    `newtxtext` + `newtxmath` (Times-compatible; `src/v2/figure_typography.py`).
    Frontiers draft uses default `article` (Computer Modern). Site type should be
    **print-adjacent serif on screen**, not a CM webfont clone (hairline at UI sizes)
    and not Inter/Roboto/system/Space Grotesk.
  - **PeerJ upload numbering ≠ R stem order:** `FIGURE_MAP` in
    `src/v2/figure_typography.py` maps `fig10_coverage_qc` → `Figure1.pdf`,
    `fig1_truth_construct` → `Figure2.pdf`, … `fig12_protocol_pass_matrix` →
    `Figure12.pdf`. `fig_study_design` is in both PDFs but **not** in `FIGURE_MAP`.
    Frontiers in-PDF order starts with study design as Figure 1 and ends coverage QC
    as Figure 13, plus appendix `fig_ext{1,2,3}_*`.
  - **Erratum (2026-08-13):** In-PDF number for `fig10_coverage_qc` is Figure 13 in **both** journal SoT manuscripts; `FIGURE_MAP` `Figure1.pdf` is a stale PeerJ upload filename, not PeerJ Figure 1.
  - **Stale pointer:** `docs/SCREG_EVAL_PROTOCOL.md` line 3 cites “Deliverable #1 of
    DESIGN.md §8”. That historical §8 is **absent** from this tree. This file is a
    GitHub Pages design contract, not a revival of that section numbering.
- Assumptions (proceed; listed as open questions where they can change):
  - English-only v1.
  - Light paper-first; no dark-mode-only-first.
  - Static GitHub Pages (HTML + CSS, optional Jekyll). No SPA unless a later
    question forces client routing.
  - Project Pages `baseurl` `/scfm-reg-audit`.
  - Site consumes **exported** PDFs/PNGs; it must not edit `paper/manuscript.tex`,
    Frontiers `manuscript.tex`, or `paper/make_figs.R`.

## Brand

- Personality: Quiet academic lab / project page. Typeset-adjacent, methods-honest,
  slightly austere. Reads like a well-set journal materials page, not a startup, not
  a lab-photography brochure, not a dark analytics console.
- Trust signals (show facts, not slogans):
  - Frozen panel: 446 TFs × 1,200 genes.
  - Dual-null (gene-label Mantel-style + degree-preserving row shuffle, \(N=999\)).
  - Protocol-pass **0/13** under the predeclared conjunction (SAP §4 / Table 4 /
    Fig. 12). Dual-null Support under full spec is **not** recovery.
  - Capsule validates published artifacts; full production rerun needs external
    data (`docs/FULL_RERUN.md`).
  - Dual license (MIT code / CC BY 4.0 manuscript, figures, derived results).
  - Zenodo concept DOI `10.5281/zenodo.21724336` and GitHub remote.
  - ORCID `0009-0001-8329-0108`; affiliation as in the author block.
- Avoid:
  - Startup landing, hero video, fake testimonials, “SOTA”, “beats”, “revolutionary”,
    purple gradients on white, Inter/Roboto, cream/terracotta hospitality editorial
    (Opus default: ~`#F4F1EA` + Fraunces/Playfair + amber — **wrong register here**).
  - Global positive or global negative claims about scFMs; “regulatory truth” /
    causal ceiling; model leaderboards; retired bootstrap CI / MDE / implied-alpha.
  - Dark-mode-first chrome, glassmorphism, auto-playing motion, stock science GIFs.
  - Collapsing PeerJ and Frontiers into one “the paper” when numbering and
    appendix set differ.

## Product goals

- Goals:
  - Give **readers** of either paper a calm place to open the matching PDF, see
    shared figures with venue-correct numbers, and copy citation metadata.
  - Give **reviewers** the protocol/SAP definitions (Support vs protocol-pass),
    claim-register pointers, and licensing without hunting the git tree.
  - Give **reproducers** the capsule validate path (`python validate_artifacts.py`)
    and an honest full-rerun boundary (this tree is not a data lake).
  - Keep **two independent paths** (PeerJ CS ↔ Frontiers Genetics) first-class.
- Non-goals:
  - Implementing the site in this design pass.
  - Editing manuscripts, `make_figs.R`, or TikZ fragments.
  - Hosting H5AD/weights/NPZ; pretending the capsule is a full rerun.
  - Interactive model explorer, live stats recalculation, or a SPA app.
  - Marketing conversion, mailing list, or “request a demo”.
  - Overwriting `.omx/plans/prd-github-pages.md` if a planner writes it.
- Success signals:
  - A reviewer can reach the correct PDF and Fig. 12 / protocol-pass table in
    two clicks from the hub, with 0/13 stated in plain language.
  - A reproducer can run the validator from the reproduce page without being
    told to download missing H5AD.
  - Venue figure numbers never silently swap (PeerJ `Figure1.pdf` ≠ Frontiers
    Figure 1).
  - No slogan appears where a gate count belongs.

## Personas and jobs

- Primary personas:
  1. **Paper reader** — has the PDF or a citation; wants figures, the other
     venue’s PDF, and a BibTeX/DOI block.
  2. **Reviewer / editor** — needs methods honesty: estimand, dual-null,
     protocol-pass conjunction, what is *not* claimed, licenses, AI-disclosure
     pointer (`paper/submission_peerj/ai_code_disclosure/`).
  3. **Capsule reproducer** — clones GitHub; wants validate + unit tests +
     FULL_RERUN pointers + accession list, not a GUI for `run_fixed_panel_audit.py`.
- User jobs:
  - “Open the PeerJ (or Frontiers) PDF.”
  - “Show me Fig. 12 / the protocol-pass matrix and what 0/13 means.”
  - “Download this figure as used in *this* venue.”
  - “Cite the capsule / preprint archive.”
  - “What can I rerun from a fresh clone?”
- Key contexts of use:
  - Desktop / laptop, journal-reviewer or literature-skimming session.
  - Possibly a phone in a corridor: PDF link and honesty strip must still work;
    figure galleries may stack.
  - Slow or blocked GitHub raw: prefer Pages-hosted small HTML/CSS; large PDFs
    may live on GitHub Releases or Zenodo (open question).

## Information architecture

Multiple paths, not a funnel. The hub is a **fork**: two venue doors plus shared
resources. Do not hide the second venue behind a dropdown-only control.

- Primary navigation (persistent, text, not icon-only):
  - Hub
  - PeerJ CS
  - Frontiers Genetics
  - Figures
  - Protocol
  - Reproduce
  - Cite & licenses
  - GitHub (external)
- Core routes/screens (pretty URLs; `.html` fallbacks if not using Jekyll):

  | Route | Purpose |
  | --- | --- |
  | `/` | Hub: one-paragraph what-this-is; honesty strip (panel, dual-null, protocol-pass 0/13); two venue cards; links to GitHub + Zenodo. No hero video, no CTA banner. |
  | `/peerj/` | PeerJ CS path: review PDF pointer, upload figure list using **`FIGURE_MAP` names** plus stable stems, tables 1–8, package README notes, AI-disclosure pointer. State that `fig_study_design` appears in the compiled PDF but is not a `FigureN.pdf` in `flat_upload/`. |
  | `/frontiers/` | Frontiers Genetics path: review PDF pointer, in-PDF Figures 1–13 (study design first, coverage QC last), Appendix A1–A3 (`fig_ext1`–`fig_ext3`), Vancouver vs PeerJ name-year noted only as a bibliographic fact — not a selling point. |
  | `/figures/` | Shared gallery indexed by **stable stem**. Each card shows PeerJ upload name (if any) and Frontiers figure/appendix number. |
  | `/figures/{stem}/` | Single-figure page: caption excerpt, venue numbers, download/view, “same asset, different number”. Stems listed below. |
  | `/protocol/` | Human summary of protocol + SAP: estimand, two nulls, protocol-pass gates, link to `docs/SCREG_EVAL_PROTOCOL.md` and `docs/SCREG_EVAL_SAP.md`. Include the Fig. 12 pass/fail legend as the visual key. |
  | `/reproduce/` | Layer 1: `validate_artifacts.py` + `unittest`. Layer 2–3: `docs/FULL_RERUN.md` is **out of capsule**. Fail-closed: missing H5AD is expected. Accessions: GSE174367, 10x PBMC multiome, GSE206767 (construct-only). |
  | `/resources/` | PDFs, JSON derivatives (authoritative `results/v2/*` names only), MANIFEST/SHA-256, environment lock pointer. No retired JSON as “results”. |
  | `/cite/` | `CITATION.cff` fields, Zenodo DOI, dual license, ORCID, correspondence email as already public in the manuscripts. |

- Content hierarchy:
  1. Which path (PeerJ vs Frontiers vs shared).
  2. What the audit **is** (fixed-panel unusualness, not ranking).
  3. Artifacts (PDF, figures, tables).
  4. How far a clone can go (validate vs full rerun).
  5. How to cite and reuse (licenses).
- Stable figure stems (canonical IDs for `/figures/{stem}/`):

  | Stem | PeerJ `FIGURE_MAP` | Frontiers in-PDF (observed order) |
  | --- | --- | --- |
  | `fig_study_design` | not in upload map (in PDF) | Figure 1 |
  | `fig1_truth_construct` | Figure2 | Figure 2 |
  | `fig2_cross_tissue_decomp` | Figure3 | Figure 3 |
  | `fig3_primary_audit` | Figure4 | Figure 4 |
  | `fig4_usability_check` | Figure5 | Figure 5 |
  | `fig5_null_diagnostics` | Figure6 | Figure 6 |
  | `fig6_spec_sensitivity` | Figure7 | Figure 7 |
  | `fig7_pertype_descriptive` | Figure8 | Figure 8 |
  | `fig8_injection_ladder` | Figure9 | Figure 9 |
  | `fig9_tf_probe` | Figure10 | Figure 10 |
  | `fig11_third_tissue_transfer` | Figure11 | Figure 11 |
  | `fig12_protocol_pass_matrix` | Figure12 | Figure 12 |
  | `fig10_coverage_qc` | **Figure1** | Figure 13 |
  | `fig_ext1_construct_mantel` | — | Appendix A1 |
  | `fig_ext2_baselines_collectri` | — | Appendix A2 |
  | `fig_ext3_honesty_policy` | — | Appendix A3 |

  Always label with the stem in `code` style on the figure page so numbering
  collisions cannot be mistaken for different science.

## Design principles

- Principle 1: **Honesty before headline.** Lead with the frozen panel, dual-null,
  and protocol-pass 0/13. Dual-null Support under the full spec is a secondary
  count with the degree-hazard sentence attached — never a victory metric.
- Principle 2: **Two papers, one audit, two doors.** PeerJ CS and Frontiers
  Genetics share stems and JSON; they do not share figure numbers or appendix
  sets. Navigation makes both paths equal.
- Principle 3: **Print-adjacent, screen-legible.** Serif body, high contrast on
  paper-white, figure aqua/grey reused as semantic pass/fail — not as a brand
  splash. UI chrome stays quiet.
- Principle 4: **Materials, not a product launch.** Every page answers “open,
  download, cite, validate.” No conversion copy. External GitHub/Zenodo links
  are utilities.
- Principle 5: **Stem is the identity; venue number is a label.** Never title a
  gallery item “Figure 1” without saying which venue. Prefer `fig10_coverage_qc`
  as the durable key.
- Tradeoffs:
  - Static Pages over SPA: slower to “app-like”, faster to trust and archive.
  - Light-only v1 over dark theme: matches PDFs; dark can wait (open question).
  - Linking large PDFs vs inlining: prefer one obvious PDF button per venue
    rather than embedding full manuscripts in HTML.
  - Showing 0/13 on the hub vs “burying the lede”: we **show** it; we do not
    editorialize it into a field-level negative (manuscript conclusions forbid that).

## Visual language

Print-adjacent academic materials. **Override** of cream/serif-italic/terracotta
editorial hospitality. Palette is taken from the papers’ own figure tokens plus
ink-on-paper, not from a magazine moodboard.

- Color (CSS variables; implementers must use these names):

  | Token | Hex | Use |
  | --- | --- | --- |
  | `--paper` | `#FAFBFC` | Page background (cool paper, **not** `#F4F1EA`) |
  | `--paper-inset` | `#FFFFFF` | Cards, figure wells |
  | `--ink` | `#1A1A1A` | Body and headings (Fig 12 header ink) |
  | `--muted` | `#5A6570` | Meta, captions, conjunction-column slate |
  | `--rule` | `#D8DCE0` | Hairline rules, table borders |
  | `--pass` | `#1B8F75` | Protocol pass / dual-null-supported fill (AQUA) |
  | `--fail` | `#E0E0E0` | Protocol fail (`grey88`) |
  | `--link` | `#1B4F72` | Text links (ink-teal, Frontiers-like `blue!50!black` spirit) |
  | `--link-hover` | `#2A78D6` | Hover only (figure BLUE; do not paint large surfaces) |
  | `--warn` | `#D99A00` | Caution callouts (degree hazard, retired artifacts) — sparse |
  | `--danger-text` | `#9B2C2C` | Error text only; do **not** reuse figure RED `#D84A4A` as decoration |
  | `--focus` | `#1B8F75` | Focus ring |

  Pass/fail chips: aqua fill + white or ink label only if contrast ≥ 4.5:1;
  otherwise aqua border + ink text on `--paper-inset`. Fail chips: grey fill is
  decorative; the word “fail” is ink, not grey-on-grey.

- Typography:
  - Headings and body: **STIX Two Text** (scientific Times successor; PeerJ
    newtx-adjacent, readable at 16–18px). Fallback: `Source Serif 4`, `Georgia`,
    `Times New Roman`, serif.
  - UI chrome (nav, badges, table headers, honesty-strip labels): **Atkinson
    Hyperlegible** (letter-disambiguation for q/ρ/TF IDs). Fallback: `Source
    Sans 3`, sans-serif. Not Inter, Roboto, Arial-as-brand, or Space Grotesk.
  - Hashes, accessions, commands: **IBM Plex Mono**. Fallback: `ui-monospace`,
    `Source Code Pro`, monospace.
  - Scale (approx): body 18px / 1.55; h1 2rem; h2 1.35rem; caption 0.9rem;
    nav 0.95rem. Max measure ~68ch for prose; figure pages may go ~80rem wide.
  - Do not load Computer Modern / Latin Modern as the UI face (too thin on
    screens). PDFs remain CM or Times as compiled.

- Spacing/layout rhythm:
  - Base 8px; section padding 32–48px; hub max-width ~72rem; figure well
    full-bleed to ~80rem.
  - Vertical rhythm like a journal: rule, small-caps or tracked label, then
    title — not a giant gradient hero.
  - Honesty strip: one horizontal band under the header (stack on small
    screens) with three facts, equal visual weight.

- Shape/radius/elevation:
  - Radius 2–4px (near-print). No 16px “app cards”.
  - Elevation: 1px `--rule` border, no drop-shadow stacks. Figure wells: inset
    border only.

- Motion:
  - Almost none. `:hover` underline or link color 150ms ease.
  - No page-load choreography, no staggered hero reveals, no count-up
    animation on “0/13” (the zero must not feel like a scoreboard).
  - `prefers-reduced-motion: reduce` → instant states.

- Imagery/iconography:
  - Primary images = paper figures (TikZ/PDF/PNG exports). Do not redraw Fig 12
    in CSS as a decorative hero.
  - No stock DNA helices, robot brains, or GPU-glow thumbnails.
  - Icons: text + simple 16px strokes if needed (PDF, GitHub). No illustrated
    icon set.
  - Optional tiny wordmark: “scReg-Eval” in STIX Two, weight 600; no logo mark
    required in v1.

## Components

- Existing components to reuse:
  - None in-repo (no frontend). Reuse **figure color tokens** and **caption
    language** from the manuscripts. Reuse `CITATION.cff` for the cite block.
- New/changed components (site-only; names are implementation hints):
  - `SiteHeader` — wordmark, primary nav, GitHub text link, venue indicator
    when on `/peerj/` or `/frontiers/`.
  - `HonestyStrip` — three facts: `446 × 1,200` panel; dual-null; protocol-pass
    `0/13`. Each fact is a label + value, not a KPI sparkline.
  - `VenueCard` — journal name, short path descriptor (“PeerJ CS report” /
    “Frontiers Genetics conversion + appendix”), PDF button, “figures in this
    numbering” link. Equal width on hub.
  - `FigureCard` / `FigurePage` — thumbnail or PDF object, stem, dual numbers,
    caption excerpt, download. Empty well if export missing.
  - `PassFailChip` — aqua/grey; always includes the word pass/fail (color is
    not the only channel).
  - `ScopeCallout` — “instance scope, not a field-level negative” (Fig 12B
    language). Warn token, not alarm red.
  - `ReproduceLayers` — numbered 1 / 2 / 3 matching `docs/FULL_RERUN.md`.
  - `CiteBlock` — preformatted citation + copy button (clipboard API; fallback
    select-all).
  - `SiteFooter` — MIT / CC BY 4.0, DOI, correspondence, “capsule v0.3.0”.
- Variants and states:
  - VenueCard: PeerJ | Frontiers (typographic difference only, same chrome).
  - PassFailChip: pass | fail | not-applicable (appendix-only stems on PeerJ
    path).
  - FigureCard: available | missing-export | PeerJ-upload-only |
    Frontiers-appendix-only.
- Token/component ownership:
  - Visual tokens: this `DESIGN.md`.
  - Numerical claims: `docs/SCREG_EVAL_SAP.md`, authoritative JSON, manuscripts.
    The site must not invent counts. If SAP and a page disagree, SAP wins and
    this file gets an open question.
  - Figure pixels: `paper/make_figs.R` / `paper/make_figs_extension.R` outputs;
    site does not restyle plot internals.

## Accessibility

- Target standard: WCAG 2.2 AA on the static pages.
- Keyboard/focus behavior:
  - All nav, PDF, and copy controls reachable in DOM order.
  - Visible focus ring `2px solid var(--focus)` with 2px offset; never
    `outline: none` without a replacement.
  - Skip link “Skip to content” as first focusable.
- Contrast/readability:
  - Ink on paper ≥ 12.5:1 body; muted on paper ≥ 4.5:1.
  - Aqua `#1B8F75` on white is in the ~4.5:1 region for large text; do not use
    aqua for small body text. Fail grey is never the only indicator.
  - Caption length and 18px body beat “tiny journal 9pt on screen”.
- Screen-reader semantics:
  - One `h1` per page; venue names in text, not color.
  - Honesty strip as a list or description list.
  - Figures: `alt` = short stem + venue number + one-clause caption; long
    caption in adjacent text (not only in `alt`).
  - Copy-citation button: live region “Copied citation.”
- Reduced motion and sensory considerations:
  - Honor `prefers-reduced-motion`.
  - No autoplay; no contrast-only pass/fail.
  - Do not rely on aqua vs grey for color-blind users — include “pass”/“fail”
    text (already in Fig 12 legend).

## Responsive behavior

- Supported breakpoints/devices:
  - `≤640px` phone; `641–960px` tablet; `≥961px` desktop. No separate app.
  - Browsers: last two Chrome/Firefox/Safari/Edge; no IE.
- Layout adaptations:
  - Nav: horizontal wrap or a details/summary “Menu” on small screens (no
    hamburger-only mystery).
  - Hub venue cards: one column → two columns.
  - Honesty strip: three columns → stacked rows.
  - Figure pages: image full width of well; dual-number meta stacks under the
    title.
  - Tables on `/protocol/`: horizontal scroll inside a well rather than
    squashing 13×7 gates to unreadability.
- Touch/hover differences:
  - Hit targets ≥ 44px for PDF and copy.
  - Hover underline is extra; links already underlined or marked “PDF”.

## Interaction states

- Loading: static site — HTML should be complete without JS. If a PDF `<object>`
  is used, show a textual “Download PDF” sibling immediately (no spinner-only).
- Empty: missing figure export → well with stem name + “not exported on this
  path” + link to GitHub blob if known. Do not show broken-image icons as
  science.
- Error: 404 page in the same chrome: “No such materials path” + hub link.
  GitHub Pages 404 is otherwise unbranded — ship `404.html` if using custom
  HTML.
- Success: copy-cite confirmation only. Do not toast “Welcome!” on load.
- Disabled: N/A controls (e.g. PeerJ download for `fig_ext1`) visible but
  disabled with reason “Frontiers appendix; not in PeerJ upload map”.
- Offline/slow network: no service worker in v1. Pages must remain readable
  if webfonts fail (system Georgia/Times already in the stack). Do not hide
  content behind font-display blocking; use `font-display: swap`.

## Content voice

- Tone: The manuscripts’ voice. Precise, hedged where the SAP hedges, short
  sentences. First person plural only if quoting the paper; the site itself
  speaks in third person (“This capsule…”, “The audit finds…”).
- Terminology (use these strings; do not “simplify” into marketing):
  - **scReg-Eval** — protocol + capsule name.
  - **fixed-panel audit** — not “benchmark leaderboard”.
  - **regulatory-potential proxy** — not “ground truth” or “causal GRN”.
  - **dual-null Support** — both named randomizations after BH in the same spec.
  - **protocol-pass** — predeclared conjunction; frozen **0/13**.
  - **capsule** vs **full rerun** — never conflate.
  - Tissues: brain (GSE174367), PBMC multiome, fibroblast mix (GSE206767,
    construct-only).
- Microcopy rules:
  - Hub subhead example: “A fixed-panel audit of regulatory alignment in scRNA
    foundation-model gene graphs. Two report PDFs; one capsule.”
  - Forbidden on the site (from `docs/PAPER_OUTLINE.md` plus this brief): SOTA,
    regulatory truth, causal ceiling, “scFMs encode regulation”, “no scFM
    encodes regulation”, MDE, implied alpha, TF-block bootstrap CI as current
    uncertainty, testimonials, “request access”.
  - When stating 7/13 dual-null under **full** spec, the same sentence (or the
    next) must mention non-degree 0/13 and/or protocol-pass 0/13.
  - Buttons: “Open PeerJ PDF”, “Open Frontiers PDF”, “Download figure”, “Copy
    citation”, “View protocol” — not “Get started” or “Learn more”.
  - Dates as ISO or “13 Aug 2026”; capsule version `v0.3.0`.

## Implementation constraints

- Framework/styling system:
  - GitHub Pages static hosting. Prefer **hand-written HTML + one CSS file**
    (custom properties from Visual language) or **Jekyll with a custom layout
    and no default Minima look**. No React/Vite/Next unless a later decision
    documents why static HTML failed.
  - Project Pages require `baseurl: /scfm-reg-audit` (or equivalent) on every
    asset URL.
  - No existing component library to extend.
- Design-token constraints:
  - Pass/fail colors must match figure AQUA / grey88.
  - Do not introduce a second green for “success”.
  - Do not restyle figure PNG/PDF internals with CSS filters.
- Performance constraints:
  - HTML/CSS first; JS only for copy-citation.
  - Webfonts: subset STIX Two + Atkinson + Plex Mono; three families max.
  - Thumbnails for gallery; full PDF/PNG on figure pages. Avoid shipping every
    TikZ PDF into the hub.
  - GitHub Pages soft file-size limits: if manuscript PDFs are large, link to
    GitHub `blob`/`raw` or Zenodo instead of duplicating (open question).
- Compatibility constraints:
  - Must not edit `paper/manuscript.tex`,
    `paper/submission_frontiers_genetics/manuscript.tex`, or `paper/make_figs.R`.
  - Must not reintroduce retired bootstrap/MDE artifacts as live results
    (`docs/reports/table-enrichment-plan.md` Forbidden list; `docs/LEGACY_INFERENCE_NOTE.md`).
  - Site tree should sit in a Pages source that `validate_artifacts.py` can
    treat as local overlay if needed (see `LOCAL_WORKTREE_PREFIXES` in
    `validate_artifacts.py`) — prefer `docs/site/` or a dedicated `pages/`
    path agreed in the PRD, not scattering HTML through `paper/`.
  - Dual license in footer. Do not relicense vendor code.
- Test/screenshot expectations:
  - Manual: hub, both venue pages, one shared figure (`fig12_protocol_pass_matrix`),
    reproduce page, 404, 320 / 768 / 1280 widths, keyboard-only pass.
  - Visual Ralph is **out of scope** until a visual reference is approved.
    This file is the design contract, not a pixel baseline.
  - No automated visual regression in-repo today.

## Open questions

- [ ] Pages source location (`docs/`, `docs/site/`, `pages/`, or `gh-pages`
      branch) / owner: implementer + repo owner / impact: clone layout and
      `validate_artifacts.py` overlay list.
- [ ] Jekyll vs plain HTML / owner: implementer / impact: permalinks and
      `baseurl` plumbing; DESIGN.md allows either if the visual language holds.
- [ ] Where large PDFs live (Pages copy vs GitHub vs Zenodo) / owner: repo
      owner / impact: reviewer one-click vs repo size.
- [ ] Dark mode later? v1 is light-only; a later opt-in must keep pass/fail
      semantics and ≥ AA contrast / owner: design refresh / impact: tokens.
- [ ] Publication status labels (draft / submitted / published) per venue /
      owner: author / impact: hub wording; do not invent “published in X”
      before it is true.
- [ ] Custom domain vs `*.github.io` / owner: author / impact: canonical URLs
      in CITATION later.
- [ ] Whether PeerJ in-PDF float order should be documented beside `FIGURE_MAP`
      upload names on `/peerj/` (study design present in PDF, absent from
      `FigureN.pdf`) / owner: implementer / impact: reviewer confusion — this
      draft already requires both labels.
- [ ] Thumbnail pipeline: commit PNG wrappers vs generate in CI from TikZ /
      owner: implementer / impact: must not require editing `make_figs.R`;
      visual-QA PNGs today are scratch, not production assets.
- [ ] Chinese mirror or bilingual UI / owner: author / impact: nav + PDFs stay
      English unless explicitly requested.
- [ ] Relationship to `.omx/plans/prd-github-pages.md` when it appears / owner:
      planner vs designer / impact: PRD may pin routes or hosting; reconcile
      without overwriting either file.
- [ ] Stale `docs/SCREG_EVAL_PROTOCOL.md` “DESIGN.md §8” citation / owner:
      docs / impact: protocol readers; out of scope to rewrite in this pass
      unless asked.
