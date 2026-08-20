# GitHub Pages source (`site/`) — Hugo project

Public routes: `/` `/figures/` `/reproducibility/`.
Production URL: `https://PeterPonyu.github.io/scfm-reg-audit/`.

Hugo extended **0.165.0** (pinned in `.github/workflows/pages.yml`). The hub is a door into the audit object: frozen 446 × 1,200 panel, dual-null Support, protocol-pass 0/13.

## Local

```bash
# needs Hugo extended (CI pins 0.165.0)
export HUGO=/path/to/hugo   # optional
bash site/test.sh
python3 -m http.server 4173 --directory site/public
```

`site/test.sh` runs `site/assemble.sh` then `site/test.py`. Assemble fail-closes if the publish dir contains PDFs. Figure previews are PNG; full manuscripts are not part of the public tree or the Pages artifact.

Do not commit `site/public/`. Fonts are self-hosted woff2 under `static/fonts/` (Geist SIL OFL; no Google Fonts CDN).

## CI

Push to `main` or `workflow_dispatch` builds with pinned Hugo extended, `hugo --minify`, runs `site/test.sh`, uploads `site/public/` via `upload-pages-artifact`, deploys with `deploy-pages`.
