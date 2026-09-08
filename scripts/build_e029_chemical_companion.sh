#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-$(pwd)}
if [[ $# -gt 0 ]]; then shift; fi
AIDD_PY=${AIDD_PY:-python}
AIDD_DB=${AIDD_DB:-/mnt/medchem_taltio/wrk/yu_agent/runtime/registry/aidd.sqlite3}
QT9_REPORT=${QT9_REPORT:-$ROOT/data/e020_qt9_validation/validation-report.json}
E029_OUTPUT_ROOT=${E029_OUTPUT_ROOT:-/mnt/local/hand/yuzhang/aidd/chemical-companion-v1}
WORKERS=${WORKERS:-16}

for path in "$AIDD_DB" "$QT9_REPORT"; do
  [[ -f "$path" ]] || { echo "Required E029 input not found: $path" >&2; exit 66; }
done

INDEX_CATALOG=$($AIDD_PY -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["result"]["index_catalog"])' \
  "$QT9_REPORT")
ARTIFACT_CATALOG=${ARTIFACT_CATALOG:-$($AIDD_PY -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["artifact_catalog"])' \
  "$INDEX_CATALOG")}

mkdir -p "$E029_OUTPUT_ROOT"
$AIDD_PY -m aidd_agent.cli build-chemical-companion \
  --db "$AIDD_DB" \
  --library LIB-AFA68EE6888C \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --output-root "$E029_OUTPUT_ROOT" \
  --workers "$WORKERS" \
  --max-moving-atoms 12 \
  "$@"

echo "Chemical companion: $E029_OUTPUT_ROOT/catalog.json"
