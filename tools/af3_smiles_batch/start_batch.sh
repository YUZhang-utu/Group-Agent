#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON:-python}"
MODE="${1:---check}"
if [[ "$MODE" != --check && "$MODE" != --run ]]; then
  echo 'Usage: bash start_batch.sh [--check|--run]' >&2
  exit 2
fi
OUTPUT="$ROOT/jobs_$(date +%Y%m%d_%H%M%S)_$$"
"$PYTHON_BIN" "$ROOT/prepare_batch.py" \
  --base "$ROOT/myprotein.json" --ligands "$ROOT/ligands.tsv" \
  --ligand-copies "${AF3_LIGAND_COPIES:-1}" --output "$OUTPUT"
"$PYTHON_BIN" "$ROOT/run_batch.py" "$OUTPUT" "$MODE"
