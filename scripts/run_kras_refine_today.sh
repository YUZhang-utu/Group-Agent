#!/usr/bin/env bash
set -euo pipefail
DIRECT="${1:-$HOME/kras-direct-count-v1}"
EXAMPLES="${2:-$HOME/kras-soft-clash-v1/report.json}"
OUTPUT="${3:-$HOME/kras-refine-v1}"
HOURS="${4:-4}"
PYTHON_BIN="${PYTHON:-python}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
if [[ ! -f "$OUTPUT/gate/report.json" ]]; then
  "$PYTHON_BIN" scripts/run_kras_contact_refine.py gate --direct "$DIRECT" \
    --examples "$EXAMPLES" --output "$OUTPUT/gate"
fi
"$PYTHON_BIN" scripts/run_kras_contact_refine.py batch --direct "$DIRECT" \
  --gate "$OUTPUT/gate" --output "$OUTPUT/batch" --hours "$HOURS" --workers 20 --count 1000000
EXPORT="$OUTPUT/export-$(date -u +%Y%m%dT%H%M%SZ)"
"$PYTHON_BIN" scripts/export_kras_soft_rank.py --ranking "$OUTPUT/batch" \
  --output "$EXPORT" --allow-budgeted --count 1000000 --batch-size 100000
echo "Export report: $EXPORT/report.json"
