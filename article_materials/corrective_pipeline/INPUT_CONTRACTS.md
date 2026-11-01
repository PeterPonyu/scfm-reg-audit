# Input contracts

Resource entries contain `path` and a lowercase SHA256 of the exact local file. Examples use semantic resource names and omit personal file paths. `local_resources.json`, true inputs, runs and approved local authorization records are ignored. Hashes certify identity, not biological validity. Inputs with altered hashes must be registered before execution and cannot silently replace a frozen input.

## Independent donor extension

`multiome` must be an H5AD with explicit integer raw counts at `layers/counts`, paired GEX/ATAC rows, nonmissing `DonorID` and `cell_type`, and `Modality=multiome`. Features must have unique symbols or peak coordinates and an explicit `feature_types` mapping. A recorded reference declaration must agree with the genome reference. Donor IDs are exact; batch, site and cell type are not substituted for donor identity.

`gene_coordinates` is a tab-separated five-column file without a header: chromosome, zero-based start, half-open end, strand, symbol. The existing annotation's first occurrence of each symbol defines the locus; promoter extension is fixed before outcomes. `genome_fasta` and `genome_fai` must be the corresponding existing reference and index. `meme` contains valid MEME-format probabilities and TF identifiers. `original_panel` contains a `genes` list for reporting overlap, not for choosing successful new genes.

Model dictionaries, config, vocabulary and weights must be complete and hash-verified. The legacy Geneformer loaders use the standard basenames in one dictionary directory and the config/weights in one checkpoint directory; scGPT config, vocabulary and weights also share their standard directory. The runner checks those actual loader locations so a declared resource cannot point elsewhere. scFoundation uses an explicit gene index. UCE requires the explicit ESM2 source, multi-species chromosome table and human species offset. Only local loading is allowed.

## Brain RNA source comparison

The preferred route supplies `source_left` and `source_right`, each with `path`, `sha256`, `matrix` and `scale`. Scale is exactly one of `counts`, `cp10k` or `log1p_cp10k`. Normalization uses each input's full RNA universe before selecting fingerprint genes. Exact cell and feature keys must have a documented correspondence; the runner never strips barcode suffixes or invents sample aliases. The comparison checks full identity subsets and a deterministic 128-cell by 512-gene numeric fingerprint. It cannot certify every matrix entry or reconstruct an absent preprocessing history.

The alternative `source_candidate` is an NPZ with `left_cells`, `left_genes`, `left_matrix`, `right_cells`, `right_genes` and `right_matrix`, all readable with `allow_pickle=False`. Matrices are already on the explicitly justified comparable scale.

Optional `source_donors` is a CSV with `cell,donor`. It must uniquely map every left-input cell, using composite sample-plus-barcode keys where needed. Duplicate, missing or ambiguous keys stop the join. `source_chain` is JSON with `accession`, `source_evidence_sha256` and a sequence of `steps`, each containing `input_sha256`, `output_sha256` and the actual `command`. Adjacent steps must join by their exact hashes. Structure checks alone do not prove commands were executed: source evidence, endpoint agreement, donor mapping and replay still require review. All source-comparison reports retain `provenance_closed=false`.

## Checkpoint training overlap

`training_inventory` is a JSON list with one row per checkpoint/model. Each row has `model`, `checkpoint_sha256`, `source_evidence_sha256`, `dataset_ids` and, only when justified, `complete_dataset_inventory=true`. `evaluation_identity` is JSON keyed by the same model names, with `checkpoint_sha256`, `dataset_ids` and an evidence-backed `identity_crosswalk_complete` flag. Checkpoint binding and complete accession alias resolution are necessary even for exclusion at dataset resolution.

An overlap can be documented when the identities and checkpoint agree. An incomplete training listing or crosswalk remains unknown. Dataset-level exclusion still does not establish donor/cell-level exclusion; the runner never emits a never-pretrained claim from a model card's silence.

## Tissue-matched external assay labels

`external_labels` is a CSV with `donor,tf,target,tissue,label_status,evidence_id`. Frozen evaluation donors and genes are matched exactly. The default tissue is `bone_marrow`. Labels must come from a separately justified assay, with statuses `positive`, `verified_negative` or `unlabeled`. At least five positives and five comparators per TF/donor are required. Duplicate keys and wrong tissue are rejected before panel filtering.

Verified-negative classification and positive-versus-unlabeled occupancy are different tasks. If a TF/donor has verified negatives, its verified classification excludes unlabeled examples. Otherwise the occupancy comparison is labelled as such. Assay evidence must justify each label; adding an evidence string is not biological validation. Model and co-expression scores are paired on exactly the same scored labels. Report AUC/AP contrasts, panel and finite-score exclusions, donor coverage and correction family. Four possible contrasts per readout remain in the registered Holm family, including untestable ones at p=1.

## Biological five-criterion reference classes

`biological_readout_classes` is a JSON list of independent readout-level reference cases. Every row requires a unique `readout_id`, `reference_class` of `validated_success` or `validated_failure`, `class_evidence_id`, `class_assigned_without_gate_outcomes=true` and exactly five boolean values in `gates`: `dual_full`, `concordance`, `non_degree_sign`, `readout_sign`, `above_baseline`.

The reference class must be assigned from independent biological evidence before gate outcomes are inspected. Edge occupancy alone is insufficient. The output estimates sensitivity and specificity for each criterion and their conjunction with exact binomial intervals. Limited sample size and correlations between reference cases remain interpretive limits. Even good discrimination of the classes would not prove biological necessity of the gates.

## Randomization stress

`randomization_input` is an NPZ with a finite square `graph`, unique `tf_indices`, finite two-dimensional `controls` and one `strata` label per gene. Controls follow the exact row-major non-self TF-target ordering. The registered mechanisms are iid edge noise, zero-mean heteroskedastic shared row/column shocks and stratum-dependent mean misspecification. The last mechanism is not a type-I null. Synthetic mechanism validation cannot prove exchangeability of arbitrary real scores. Degenerate draws stop rather than being silently replaced.

## Execution acknowledgement

Commands that download data or run inference read a local run record (`authorization.json`; schema in `authorization.example.json`). It sets `approved` to true and `approved_by` to `user`, gives a free-text `instruction` describing the run, carries the SHA-256 of the current plan in `plan_sha256`, and lists the permitted items separately under `actions.experiment` and `actions.download`. Downloads also need per-resource integer `download_limits` in bytes, never above the registered limit or 5 GiB. A record that permits preparation does not permit inference, reanalysis or download. No command in this package writes the record. Downloads use the Python standard library, through the proxy named in `SC_FETCH_PROXY` or the usual proxy environment variables.
