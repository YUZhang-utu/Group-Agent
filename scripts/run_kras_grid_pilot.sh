#!/usr/bin/env bash
set -euo pipefail
PLAN="${1:-$HOME/kras-retained=v1/plan}"
OUTPUT="${2:-$HOME/kras-grid-pilot-v1}"
SAMPLE="${3:-10000}"
PYTHON_BIN="${PYTHON:-python}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
"$PYTHON_BIN" scripts/run_kras_grid_pilot.py --plan "$PLAN" --output "$OUTPUT" \
  --sample "$SAMPLE" --workers 20
echo "Grid pilot report: $OUTPUT/report.json"
