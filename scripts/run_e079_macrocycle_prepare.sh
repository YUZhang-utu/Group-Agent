#!/usr/bin/env bash
# Run from the repository root after delivering the E075-E079 working tree.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo 'Usage: bash scripts/run_e079_macrocycle_prepare.sh /absolute/new-output-root' >&2
  exit 2
fi
SOURCE=/mnt/local/hand/yuzhang/aidd/mc_data
PYTHON=${AIDD_PYTHON:-python}
OUTPUT=$1

# Validate placement before any output is created. Each stage requires fresh paths.
"$PYTHON" -c 'import sys; from pathlib import Path; s=Path(sys.argv[1]).resolve(); o=Path(sys.argv[2]).resolve(); assert s.is_dir(), "Source directory missing"; assert not o.exists(), "Output root must be new"; assert o != s and s not in o.parents, "Output must be outside source tree"' "$SOURCE" "$OUTPUT"
"$PYTHON" -m aidd_agent.macrocycle_preflight --source-dir "$SOURCE" --output "$OUTPUT/inventory"
"$PYTHON" -c 'import json,sys; from pathlib import Path; r=json.loads((Path(sys.argv[1])/"inventory/report.json").read_text()); assert not r["invalid_csv_headers"], "Review invalid CSV headers before audit"' "$OUTPUT"
"$PYTHON" -m aidd_agent.macrocycle_full_audit --source-dir "$SOURCE" --output "$OUTPUT/audit" --workers 4

# Fail closed if discovery or source bytes changed between inventory and audit.
"$PYTHON" - "$OUTPUT" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
inventory = json.loads((root / 'inventory/report.json').read_text())
audit = json.loads((root / 'audit/report.json').read_text())
expected = {row['path']: row['sha256'] for row in inventory['files']}
observed = {row['path']: row['sha256'] for row in audit['csv_sources'] + audit['mol2_sources']}
if audit['status'] != 'complete' or observed != expected:
    raise SystemExit('Source coverage/hash mismatch; review before block preparation')
PY
"$PYTHON" -m aidd_agent.macrocycle_preflight --audit "$OUTPUT/audit" --output "$OUTPUT/readiness"
