#!/usr/bin/env bash
set -euo pipefail
# Run after the operator's normal module load. This never runs a docking job.
if [[ $# -gt 1 ]]; then printf 'Usage: bash scripts/probe_docking_tools.sh [module-name]\n' >&2; exit 2; fi
if [[ $# -eq 1 ]]; then
  if ! type module >/dev/null 2>&1; then printf 'The module command is unavailable in this shell. Load your site environment first.\n' >&2; exit 2; fi
  module load "$1"
fi
"${AIDD_PY:-python}" - <<'PY'
import json, os, shutil
from pathlib import Path
root = os.environ.get('SCHRODINGER')
commands = {}
for name in ('maestro', 'glide', 'ligprep', 'prepwizard', 'structconvert', 'PLANTS', 'plants'):
    found = shutil.which(name)
    candidates = [Path(root)/name, Path(root)/'utilities'/name] if root else []
    commands[name] = found or next((str(p) for p in candidates if p.is_file()), None)
print(json.dumps(dict(schrodinger_root=root, commands=commands,
    scope='Read-only installation discovery; licenses, grids, docking execution and quality not tested'), indent=2))
PY
