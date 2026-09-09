#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-$(pwd)}
if [[ $# -gt 0 ]]; then shift; fi
AIDD_PY=${AIDD_PY:-python}
QT9_REPORT=${QT9_REPORT:-$ROOT/data/e020_qt9_validation/validation-report.json}
E029_OUTPUT_ROOT=${E029_OUTPUT_ROOT:-/mnt/local/hand/yuzhang/aidd/chemical-companion-v1}
E029_SOURCE_ROOT=${E029_SOURCE_ROOT:-/mnt/local/hand/yuzhang/aidd/mc_data}
WORKERS=${WORKERS:-16}

for path in "$QT9_REPORT"; do
  [[ -f "$path" ]] || { echo "Required E029 input not found: $path" >&2; exit 66; }
done

INDEX_CATALOG=$($AIDD_PY -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["result"]["index_catalog"])' \
  "$QT9_REPORT")
ARTIFACT_CATALOG=${ARTIFACT_CATALOG:-$($AIDD_PY -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["artifact_catalog"])' \
  "$INDEX_CATALOG")}

SOURCE_ARGS=()
EXPLICIT_SOURCE=0
for argument in "$@"; do
  if [[ "$argument" == "--source" || "$argument" == --source=* ]]; then
    EXPLICIT_SOURCE=1
  fi
done
if [[ "$EXPLICIT_SOURCE" -eq 0 ]]; then
  for name in split_0001 split_0002; do
    path="$E029_SOURCE_ROOT/$name.mol2"
    [[ -f "$path" ]] || { echo "Required E029 MOL2 source not found: $path" >&2; exit 66; }
    SOURCE_ARGS+=(--source "$name=$path")
  done
fi

mkdir -p "$E029_OUTPUT_ROOT"
$AIDD_PY -m aidd_agent.cli build-chemical-companion \
  --library LIB-AFA68EE6888C \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --output-root "$E029_OUTPUT_ROOT" \
  --workers "$WORKERS" \
  --max-moving-atoms 12 \
  "${SOURCE_ARGS[@]}" \
  "$@"

echo "Chemical companion: $E029_OUTPUT_ROOT/catalog.json"
