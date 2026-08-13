# GitHub Pages source (`site/`) — Hugo project

Public routes: `/` `/peerj/` `/frontiers/` `/figures/` `/reproducibility/`.
Production URL: `https://PeterPonyu.github.io/scfm-reg-audit/`.

Hugo extended **0.165.0** (pinned in `.github/workflows/pages.yml`). The hub is a door into the audit object and a door to materials (two PDFs, figure catalog, capsule). Venue names are file labels.

## Local

```bash
# needs Hugo extended (CI pins 0.165.0)
export HUGO=/path/to/hugo   # optional
bash site/test.sh
python3 -m http.server 4173 --directory site/public
```

`site/test.sh` runs `site/assemble.sh` then `site/test.py`. Assemble fail-closes unless both SoT PDFs exist:

- `products/peerj-manuscript.pdf` ← `paper/manuscript.pdf`
- `products/frontiers-manuscript.pdf` ← `paper/submission_frontiers_genetics/manuscript.pdf`

Do not commit `site/public/` or the product PDFs. Fonts are self-hosted woff2 under `static/fonts/` (Geist SIL OFL; no Google Fonts CDN).

## CI

Push to `main` or `workflow_dispatch` builds with pinned Hugo extended, `hugo --minify`, copies the two PDFs, runs `site/test.sh`, uploads `site/public/` via `upload-pages-artifact`, deploys with `deploy-pages`.
