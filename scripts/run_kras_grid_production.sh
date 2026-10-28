#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUTPUT="${2:-$HOME/kras-grid-production-v1}"
# OS lock is released on exit, including an interrupted job.
exec 9>"${OUTPUT}.run.lock"
flock -n 9 || { echo 'Another production run holds this output lock.' >&2; exit 1; }
"${PYTHON:-python}" scripts/run_kras_grid_production.py \
  --plan "${1:-$HOME/kras-retained=v1/plan}" \
  --output "$OUTPUT" \
  --workers "${WORKERS:-20}" --cap "${CAP_PER_TARGET:-1000000}" --batch-size 100000
