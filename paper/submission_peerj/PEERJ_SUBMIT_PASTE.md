# PeerJ Computer Science — paste into the online boxes

Copy one block at a time. Do not paste local disk paths. Suggested editors/reviewers remain OPEN (you fill those).

---

## Title

```
scReg-Eval: a fixed-panel audit of regulatory alignment in single-cell RNA foundation-model gene graphs
```

## Short title

```
Fixed-panel audit of scFM gene graphs
```

## Abstract

Paste the compiled PDF abstract (Background / Methods / Results / Conclusions). Source: `paper/manuscript.tex` `\begin{abstract}`. Do not shorten the Results sentence that jointly states dual-null Support, non-degree 0/13, and protocol-pass 0/13.

## Keywords

```
single-cell RNA-seq; foundation models; gene regulatory networks; co-expression confounding; randomization inference; benchmark audit
```

## Subject areas

```
Bioinformatics; Computational Biology; Artificial Intelligence; Data Science
```

## Article type / scope (if asked)

```
Research article. The study audits whether single-cell RNA foundation-model gene graphs align with a sequence- and accessibility-derived regulatory-potential proxy after controlling co-expression and structural covariates. It reports finite-panel randomization inference, descriptive cell-type analyses, and a pipeline-sensitivity diagnostic. The proxy is not causal regulatory ground truth.
```

---

## Author (single author)

```
Zeyu Fu
```

## Affiliation

```
State Key Laboratory of Trauma and Chemical Poisoning, Institute of Combined Injury, Chongqing Engineering Research Center for Nanomedicine, College of Preventive Medicine, Army Medical University, Chongqing, China
```

## Corresponding author

```
Zeyu Fu
fuzeyu99@126.com
```

## ORCID

```
0009-0001-8329-0108
```

## Sole-author / contributions statement

```
This is a sole-author article. Zeyu Fu conceived the study, designed the protocol, wrote the software, performed the analyses, and wrote the manuscript.

CRediT: Conceptualization; Methodology; Software; Validation; Formal analysis; Investigation; Data curation; Visualization; Writing — original draft; Writing — review and editing; Project administration.
```

## Funding

```
This work received no dedicated or external funding.
```

## Competing interests

```
The author declares no competing interests.
```

## Ethics / human-subjects

```
This study reuses public, de-identified datasets and performs no new human-subject recruitment or intervention.
```

## Data and code availability

```
Code, fixed-panel manifests, seed metadata, model/readout coverage tables, authoritative audit outputs, and figure sources are available at https://github.com/PeterPonyu/scfm-reg-audit and archived at https://doi.org/10.5281/zenodo.21724336. Public inputs: GSE174367 (brain snATAC-seq), 10x Genomics 10k PBMC multiome, and GSE206767 (fibroblast/induced-neuron ATAC pool). Large upstream datasets and model weights are not redistributed; the archive records their source identifiers and the exact manifests used.
```

## Licenses (if asked)

```
Original code: MIT. Manuscript, figures, tables, and derived public results: CC BY 4.0.
```

---

## AI use (manuscript declaration)

```
A large language model was used as an auxiliary aid for limited engineering edits to analysis code. No AI system generated scientific data, randomization outputs, or any reported number. The author takes full responsibility for the manuscript and accompanying software. No AI tool is listed as an author. Original code, the prompts used, and the resulting code are provided as Supplemental Data S1–S3.
```

## AI in computer code (PeerJ CS form)

Upload **three separate** compressed files. Do not upload the `ai_code_disclosure/` working folder.

**S1 — original code**  
File: `Supplemental_Data_S1_original_code.zip`  
SHA-256: `dea581473a6c38bf80a0da741f44dc970cecb6beb5dcc9dde376b418d6acae92`

**S2 — prompts**  
File: `Supplemental_Data_S2_prompts.zip`  
SHA-256: `5e91d71bb7002a8c4ba96e081ce6490f589474303bc1ef52b50bd982d7da050b`

**S3 — resulting code**  
File: `Supplemental_Data_S3_resulting_code.zip`  
SHA-256: `ceadeaa3a8bf2584cedc75b9ead0bd9b67b89783be45c7b235e4790eb9689af1`

Form sentence:

```
A large language model was used as an auxiliary aid for limited engineering edits to analysis code. No AI generated scientific data or reported numbers. The author takes full responsibility. Supplemental Data S1 (original code), S2 (prompts), S3 (resulting code).
```

On disk these zips are in `paper/submission_peerj/supplemental/`. Paste only the filenames into the journal form, not that path.

---

## Files to upload

Review PDF (latest PeerJ CS version; do not upload `paper/manuscript.pdf`):

```
paper/submission_peerj/flat_upload/manuscript.pdf
```

Main figures (printed 1–13). Figure 1 = study design; Figure 13 = coverage QC:

```
Figure1.pdf
Figure2.pdf
Figure3.pdf
Figure4.pdf
Figure5.pdf
Figure6.pdf
Figure7.pdf
Figure8.pdf
Figure9.pdf
Figure10.pdf
Figure11.pdf
Figure12.pdf
Figure13.pdf
```

Appendix figures:

```
FigureA1.pdf
FigureA2.pdf
FigureA3.pdf
```

Those PDFs live in `paper/submission_peerj/flat_upload/`. Optional convenience zip of the flat set (does **not** include S1–S3): `paper/submission_peerj/upload.zip`.

Checksums: `paper/submission_peerj/SHA256SUMS.txt` (flat upload) and `paper/submission_peerj/supplemental/SHA256SUMS_AI.txt` (S1–S3).

---

## Cover letter

PeerJ does not require one. If a box appears, paste:

```
Please consider this sole-author Research article for PeerJ Computer Science. scReg-Eval is a fixed-panel audit protocol and software capsule for single-cell RNA foundation-model gene graphs. Dual-null Support is not treated as regulatory recovery: protocol-pass is 0/13 under the predeclared gates. Code and data are at https://github.com/PeterPonyu/scfm-reg-audit and https://doi.org/10.5281/zenodo.21724336.
```

## Suggested editors / reviewers

Leave empty until you confirm names and conflicts (`HUMAN_GATES.md` OPEN item).
