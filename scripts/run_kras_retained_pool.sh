#!/usr/bin/env bash
set -euo pipefail
DIRECT="${1:-$HOME/kras-direct-count-v1}"
RANKING="${2:-$HOME/kras-refine-v1/batch}"
OUTPUT="${3:-$HOME/kras-retained-v1}"
PYTHON_BIN="${PYTHON:-python}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
if [[ ! -f "$OUTPUT/plan/report.json" ]]; then
  "$PYTHON_BIN" scripts/export_kras_retained_pool.py prepare \
    --direct "$DIRECT" --ranking "$RANKING" --output "$OUTPUT/plan"
fi
"$PYTHON_BIN" scripts/export_kras_retained_pool.py export \
  --plan "$OUTPUT/plan" --output "$OUTPUT/export" --batch-size 100000
echo "Retained pool export report: $OUTPUT/export/report.json"
