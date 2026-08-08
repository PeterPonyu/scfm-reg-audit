# manuscript.docx Figure Pipeline

## Current State (Aug 7, 2026)

**manuscript.docx is STALE** — contains SVG figures from Aug 3, before your refinements.

## Build Chain

```
1. paper/make_figs.R (YOUR REFINED CODE ✅ Aug 7 23:36)
   ↓ uses tikzDevice
   ↓ generates
   
2. paper/figs/fig*.tex (LaTeX/TikZ fragments ✅ Aug 7 23:36)
   ↓ compiled by
   
3. src/v2/build_peerj_package.py
   ↓ creates standalone PDFs
   ↓ 
   
4. submission_peerj/flat_upload/Figure*.pdf (✅ Ready for rebuild)
   ↓ (manual conversion?)
   ↓
   
5. submission_peerj/manuscript.docx (❌ STALE, Aug 3)
   └─ Contains: 11 embedded SVG files (rId17.svg through rId67.svg)
```

## Key Facts

1. **make_figs.R uses tikzDevice** to generate `.tex` files (LaTeX/TikZ code)
2. **Your refinements ARE in the .tex files** (✅ committed)
3. **docx contains pre-rendered SVGs** from Aug 3 (before refinements)
4. **PDF pipeline is current** — `flat_upload/manuscript.pdf` uses Aug 7 figures

## What You Need to Do

### Option 1: Rebuild Submission PDFs (Recommended)
```bash
cd paper
python3 ../src/v2/build_peerj_package.py
```

This regenerates:
- `submission_peerj/flat_upload/Figure1.pdf` through `Figure11.pdf` ← **WITH YOUR REFINEMENTS**
- `submission_peerj/flat_upload/manuscript.pdf` ← Main submission document

**PeerJ likely uses these PDFs**, not the docx.

### Option 2: Update manuscript.docx (If Required)

The docx appears to be a **secondary format**. README emphasizes PDF:
> "The current human-review PDF is `flat_upload/manuscript.pdf`"

If docx is needed:

1. **Wait for PDF rebuild** (Option 1 above)
2. **Manual conversion** (docx SVGs embedded manually, not auto-generated):
   ```bash
   # Convert PDFs to SVG
   cd submission_peerj/flat_upload
   for i in {1..11}; do
     pdf2svg Figure$i.pdf Figure$i.svg
   done
   
   # Manually re-embed in Word (or use pandoc/LibreOffice)
   ```

3. **Or skip docx entirely** if PeerJ submission only needs PDFs

## Your Refinements Status

### ✅ Successfully Applied to Source Code:
- Fig 4C: grey20 contrast fix
- Fig 11D: expanded x-axis [0.40, 0.60]
- Fig 7A/B: bootstrap confidence intervals
- Fig 6B: continuous gradient fill
- Fig 3B: closest-q annotations
- Fig 9B: construction covariate annotation
- Fig 8D: baseline reference bars
- Fig 10B: coverage percentages
- Typography: hyphen consistency, no comma separators

### ✅ Regenerated Files:
- `paper/figs/fig*.tex` (Aug 7 23:36)
- All 11 LaTeX fragments up-to-date

### ⏳ Pending Rebuild:
- `submission_peerj/flat_upload/Figure*.pdf` ← **Run build_peerj_package.py**
- `submission_peerj/flat_upload/manuscript.pdf` ← Main deliverable

### ❌ Stale (Likely Unused):
- `submission_peerj/manuscript.docx` (Aug 3, pre-refinement SVGs)

## Recommendation

**Run the PDF rebuild from the paper directory:**
```bash
cd projects/scfm-reg-audit/paper
python3 ../src/v2/build_peerj_package.py
```

This updates the actual submission files PeerJ will use. The docx is likely auxiliary/archive only.
