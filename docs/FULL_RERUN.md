# Full-rerun recipe — scReg-Eval v2

This document is the **development / experiment-workspace** full-pipeline recipe.
It is **not** something a fresh public-capsule clone can run end-to-end without
external inputs.

## 0. Two trees (do not conflate)

| Tree | Role | Full rerun? |
|------|------|-------------|
| **Publish / audit capsule** `~/Desktop/scfm-reg-audit` (GitHub `PeterPonyu/scfm-reg-audit`) | Manuscript SoT, public JSON, validator, PeerJ package | **No.** Validate only. |
| **Experiment workspace** `singlecell-genomics-research/projects/scfm-reg-audit` | Heavy data, NPZ caches, model weights, vendor code | **Yes** (this recipe). |

### Capsule fail-closed gate (fresh clone)

```bash
cd /path/to/scfm-reg-audit   # public capsule checkout
python validate_artifacts.py
python -m unittest discover -s src/tests
```

Expected: **PASS** without downloading H5AD/NPZ/weights.

The following will **fail closed** on a capsule-only tree (missing inputs, not silent
wrong results):

```bash
# These paths / artifacts are NOT redistributed in the capsule:
python src/v2/run_fixed_panel_audit.py
# → requires results/v2/*.npz graph caches and raw ATAC/RNA under SCREG_DATA_ROOT
python src/v2/build_atac_graph_v2.py
# → requires genome FASTA, ATAC H5AD, JASPAR, monorepo-scale disk
```

If `results/v2/*.npz` or `SCREG_DATA_ROOT` / `SCFM_BRAIN_ATAC` are unset, drivers
must error with a clear missing-input message. Do not invent empty graphs.

`src/v2/` in the **capsule working tree** may exist as a local development overlay
(see `validate_artifacts.py` `LOCAL_WORKTREE_PREFIXES`). The **closed release
tarball** ships only the sanitized `src/` statistical surface listed in
`MANIFEST.json`, not a pretend full monorepo.

Pinned panel hashes (protocol freeze, in-repo — not OSF):

- `manifest_sha256` = `6b203fcfab45dc600f84d2149c7f5f94e1a876f584529a0b465694e170b4f848`
- `tf_panel_sha256` = `b07ae73888cd2e075cd1992f73b5fac9a2fa9c5d8f73c596eac653d259eae8da`

## 1. Three reproducibility layers

1. **Artifact validation** — no downloads. Capsule: `python validate_artifacts.py`
   plus `python -m unittest discover -s src/tests`.
2. **Statistical reruns** — needs datasets below plus cached FM graphs in
   `results/v2/*.npz` (not redistributed; regenerate per §4 or treat as inputs).
3. **Full pipeline** — datasets + model weights + vendor code. **Experiment workspace
   only.** Commands below assume monorepo root
   `projects/scfm-reg-audit` and `ENVIRONMENT.example`.

## 2. Environment (locked 2026-07-31)

- Python 3.13.5; numpy 2.2.6; scipy 1.16.3; anndata 0.12.10; torch 2.12.0+cu130;
  scikit-learn 1.8.0; pyfaidx 0.9.0.4; pandas 2.3.3
- R 4.3.3; ggplot2 4.0.3; dplyr 1.2.1; patchwork 1.3.2; jsonlite 2.0.0; tikzDevice 0.12.6
- LaTeX: pdfLaTeX with tikz, booktabs, natbib, mdframed, placeins, helvet

## 3. Datasets (SHA pins)

| role | accession / source | local file | SHA-256 | bytes |
|---|---|---|---|---|
| brain snATAC | GSE174367 (GEO) | `GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad` | `ca6dac8097bba732f3453e9e78139ec1b11d6f63ae33e3f44d243f2377eb3ae8` | 1,739,124,186 |
| brain RNA (cross-study) | preprocessed AD cohort (project-internal preprocessing of public AD snRNA) | `ad_hm_prepped.h5ad` | `d6e8f4f2b63e13e29dc83975bd4a5789a9f742a227ccc34cb983628001df4fa0` | 80,638,259 |
| PBMC multiome RNA | 10x Genomics 10k PBMC Multiome | `pbmc10k_rna.h5ad` | `d9be197a3a777b0f97b1a2385dd18ed15d3aba9f009e36ccbef1d68261efd31d` | 128,244,328 |
| PBMC multiome ATAC | 10x Genomics 10k PBMC Multiome | `pbmc10k_atac.h5ad` | `d2afd953ac98f7626fb501024e5eab76a38be95751e62a488e06d2cfbe372a15` | 789,198,792 |
| PBMC multiome raw bundle | 10x Genomics 10k PBMC Multiome | `pbmc10k_multiome_filtered.h5` | `3897be5c916a66def9273049f6c5418ae042111f1b64ce769198253f56bf356d` | 166,323,468 |
| fibroblast reprogramming ATAC | GSE206767 (GEO; day-0 BJ + day-3 FiN pooled) | `GSE206767_filtered_peak_bc_matrix.h5ad` | `1714c1b63a66f1a121c88e46612c26819a13b0445f74a082b7c1506257af99a2` | 700,427,288 |
| JASPAR motifs | JASPAR2024 CORE vertebrates | (site archive) | — | — |
| gene coordinates | hg38 gene annotation (project-derived `gene_coords_hg38.tsv`) | data/annotation | — | — |

Accession URLs: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367`,
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206767`, 10x Genomics
datasets portal ("10k PBMC from a healthy donor, multiome"). Download dates: 2026-07.

## 4. Model weights (SHA pins)

| model | source | checkpoint | SHA-256 | license note |
|---|---|---|---|---|
| Geneformer V2-104M | HuggingFace `ctheodoris/Geneformer-V2-104M` | `model.safetensors` | `fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f` | Apache-2.0 |
| Geneformer dicts | same repo (`geneformer/` pkl dicts) | token/median/name-id dicts | `67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f` / `a51c53f6a771d64508dfaf61529df70e394c53bd20856926117ae5d641a24bf5` / `fabfa0c2f49c598c59ae432a32c3499a5908c033756c663b5e0cddf58deea8e1` | Apache-2.0 |
| scGPT-human | scGPT authors (`bowanglab/scGPT`, whole-human checkpoint) | `best_model.pt` | `6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922` | MIT |
| scGPT vocab | same | `vocab.json` | `acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f` | MIT |
| scFoundation | BioMap (`biomap-research/scFoundation`) encoder checkpoint | `model.pt` | `2446c2fb99a0e183cdc9872b7e9e55ac74f7ec82db7ba0799c83be75cde66dc6` | see upstream |
| UCE 4-layer | Stanford SNAP (`snap-stanford/UCE`) | `4layer_model.torch` | `acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560` | BSD-3 |
| UCE protein emb. | same | `human_esm2.pt` | `a210e1cc7901513999b2bca3836ba9e2f203cd008be4e9a9d6412a2267de9748` | BSD-3 |

The UCE checkpoint and ESM2 hashes are enforced at load time in
`src/v2/pbmc_uce_eval_v2.py` and `src/v2/fm_readout_uce.py`
(`EXPECTED_CHECKPOINT_SHA256`, `EXPECTED_ESM2_SHA256`). Geneformer pickle
dictionaries and the scGPT `best_model.pt` are SHA-pinned before unpickle /
`torch.load` in `src/v2/fm_readout.py`. Remaining small pickles such as
`data/uce/species_offsets.pkl` are **trusted FULL_RERUN inputs only** (not
capsule content): obtain them from the pinned upstream UCE release and do not
load untrusted pickles on the public audit path.

## 5. Ordered commands (experiment workspace)

Set `SCREG_DATA_ROOT` (datasets + models) and `SCFM_BRAIN_ATAC` first.
Working directory: monorepo `projects/scfm-reg-audit` (or an equivalent full tree
that actually contains `results/v2` NPZ caches and vendor code).

```bash
# 1. Panel + proxy (CPU, hours for motif scan)
python src/v2/freeze_gene_manifest.py
python src/v2/build_atac_graph_v2.py          # per tissue (brain, PBMC, GSE206767)

# 2. FM graphs (model-specific; GPU recommended)
python src/v2/crossmodal_v2.py
python src/v2/readout_attention_v2.py
python src/v2/insilico_ko_v2.py
python src/v2/crossmodal_scf_v2.py
python src/v2/crossmodal_uce_v2.py
python src/v2/gen_floor_attn_graphs_v2.py
python src/v2/pbmc_eval_v2.py
python src/v2/generate_pbmc_scgpt_graph_v2.py
python src/v2/pbmc_uce_eval_v2.py

# 3. Authoritative statistics (CPU)
python src/v2/run_fixed_panel_audit.py
python src/v2/brain_coexp_baseline_null.py
python src/v2/pbmc_coexp_baseline_null.py
python src/v2/subdivide_injection.py

# 4. Probe (CPU, ~10 min)
python src/v2/tf_disjoint_split.py
python src/v2/build_pair_features.py
python src/v2/run_pair_probe.py
python src/v2/pair_probe_stats.py
python src/v2/pair_probe_sensitivity.py

# 5. Figures + manuscript (publish capsule paper/ is SoT for PeerJ)
cd paper && Rscript make_figs.R && latexmk -pdf manuscript.tex

# 6. Packages (from publish capsule checkout)
python src/v2/build_peerj_package.py
python src/v2/build_release_capsule.py
```

Verification: compare regenerated JSONs against `results/v2/*.json` (field-level),
export scrubbed `*.public.json` into the capsule, then run
`python validate_artifacts.py` inside the capsule. The
`ORIGINAL_TO_PUBLIC_HASH_BRIDGE.json` records expected release hashes.

## 6. Resources

- Disk: ~6 GB datasets + ~6 GB model weights + ~2 GB caches.
- RAM: motif scan ~8 GB; FM inference 8–24 GB (GPU optional; production numbers
  in this release were produced on a 24 GB laptop GPU and CPU).
- Time: proxy construction dominates (hours); statistics ~30 min; probe ~10 min.

## 7. What is NOT redistributed (capsule)

Cached NPZ graphs (`results/v2/*.npz`), vendored model code (`src/v2/scf_vendor`,
`src/v2/uce_vendor`), model weights, raw H5AD/H5/FASTA, and retired legacy JSONs.
The public capsule validates published artifacts; this recipe regenerates them in the
experiment workspace only.
