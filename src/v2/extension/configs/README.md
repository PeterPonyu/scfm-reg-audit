# Extension configurations

These configurations describe additional tissues and graph methods. They keep the original
fixed-panel analysis separate. The original `brain` RNA label was later
corrected to ascending aorta, as described in the repository README. The loader reads JSON; the YAML files are optional copies.

| File | Contents |
|---|---|
| `tissues.json` | Tissue identifiers, modalities, analysis roles, panel rules and input paths |
| `methods.json` | Graph and method identifiers, BH-family membership and estimand descriptions |
| `*.yaml` | Optional editable copies; the loader uses JSON |

The original audit validators require the historical `brain` and `pbmc` contexts and thirteen
full-specification rows. Additional outputs go to `results/v2/extension/`; summary tables go
to `docs/reports/extension-claim-pack/`.

Entries tagged `out_of_scope`, including the Cancer/Dev RNA collections, are excluded from
Support and `G_ATAC` calculations. The historical BMMC entry remains `construct_candidate`.
Changing it to `primary_audit` requires a panel rule and a new statistical analysis plan.
This configuration is separate from the registered marrow analysis in `article_materials/`.

The construct command's `--execute` option runs Mantel and decomposition calculations on existing
`G_ATAC` NPZ files. It does not rebuild the graphs.

```bash
python src/v2/extension/cli.py registry --json
python src/v2/extension/cli.py claim-pack
```
