#!/usr/bin/env bash
set -euo pipefail
DIRECT="${1:-$HOME/kras-direct-count-v1}"
OUTPUT="${2:-$HOME/kras-soft-v1}"
PYTHON_BIN="${PYTHON:-python}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
if [[ ! -f "$OUTPUT/pilot/report.json" ]]; then
  echo "Run the 10000-conformer pilot first; see to_human/E118_KRAS_MILLION_DOCKING.md" >&2
  exit 2
fi
"$PYTHON_BIN" scripts/kras_soft_rank.py --direct "$DIRECT" --output "$OUTPUT/full" \
  --sample 0 --workers 20 --pilot "$OUTPUT/pilot/report.json"
"$PYTHON_BIN" scripts/export_kras_soft_rank.py --ranking "$OUTPUT/full" \
  --output "$OUTPUT/export" --count 1000000 --batch-size 100000
