#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON:-python}"
MODE="${1:---check}"
if [[ "$MODE" == --resume ]]; then
  [[ $# -eq 2 && -f "$2" ]] || { echo 'Usage: bash start_batch.sh --resume /path/to/jobs/run/report.json' >&2; exit 2; }
  REPORT_DIR="$(cd -- "$(dirname -- "$2")" && pwd)"
  REPORT="$REPORT_DIR/$(basename -- "$2")"
  BUNDLE="$(cd -- "$REPORT_DIR/.." && pwd)"
  exec "$PYTHON_BIN" "$ROOT/run_batch.py" "$BUNDLE" --run --resume "$REPORT"
fi
if [[ "$MODE" != --check && "$MODE" != --run ]]; then
  echo 'Usage: bash start_batch.sh [--check|--run|--resume REPORT]' >&2
  exit 2
fi
OUTPUT="$ROOT/jobs_$(date +%Y%m%d_%H%M%S)_$$"
"$PYTHON_BIN" "$ROOT/prepare_batch.py" \
  --base "$ROOT/myprotein.json" --ligands "$ROOT/ligands.tsv" \
  --ligand-copies "${AF3_LIGAND_COPIES:-1}" --output "$OUTPUT"
"$PYTHON_BIN" "$ROOT/run_batch.py" "$OUTPUT" "$MODE"
