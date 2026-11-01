"""Input locations for the Revision 3 analyses.

Each location can be set with an environment variable (SCREG_<NAME>, for example
SCREG_CORRECTIVE_RUN) or in paths.local.json next to this file (schema in paths.example.json).
Otherwise it follows the layout in which the analyses were run, or the layout of the public
deposit when this file sits in it. External data files have no default and must be set.
"""
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
W = HERE.parent
_LOCAL = HERE / 'paths.local.json'
_CFG = json.loads(_LOCAL.read_text()) if _LOCAL.exists() else {}
_DEPOSIT = (W.parent / 'corrective_pipeline' / 'screg_validation').is_dir()


def _given(name):
    value = os.environ.get('SCREG_' + name.upper()) or _CFG.get(name)
    return Path(value).expanduser() if value else None


def where(name, default):
    return _given(name) or Path(default)


def need(name):
    path = _given(name)
    if path is None:
        raise SystemExit(f'set "{name}" in {_LOCAL} or the environment variable SCREG_{name.upper()}')
    return path


REV = where('revisions', W.parent)
EXT = where('evidence_workspace', REV / 'earlier_version')
VAL = where('validation', W.parent / 'corrective_pipeline' if _DEPOSIT else EXT / 'validation')
CORR = where('corrective_run', VAL / 'corrective_run')
HIST_INPUTS = where('historical_inputs', EXT / 'build/scientific_extension/inputs')
HIST_STRESS = where('historical_stress', EXT / 'build/scientific_extension/stress')
