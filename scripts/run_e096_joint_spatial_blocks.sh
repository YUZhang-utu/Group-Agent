#!/usr/bin/env bash
set -euo pipefail
# Run from Group-Agent in the same RDKit environment as the original descriptor build.
ROOT="${1:?Supply completed final-blocks-e094 directory}"
BUILD="${2:?Supply original blocks-20000/backbone descriptor directory}"
OUT="${3:?Supply new final-blocks-e096-joint directory}"
PYTHON="${AIDD_PYTHON:-python}"
WORKERS="${AIDD_SPATIAL_WORKERS:-8}"
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$OUT/logs" "$OUT/sqlite-tmp"
export SQLITE_TMPDIR="$OUT/sqlite-tmp"
# flock releases automatically on exit or interruption; never allow two writers.
exec 9>"$OUT/run.lock"
flock -n 9 || { echo "Another E096 process holds the output lock" >&2; exit 1; }
"$PYTHON" -u -m aidd_agent.joint_spatial_profiles \
  --profiles "$ROOT/profiles" --build "$BUILD" --output "$OUT/profiles" \
  --workers "$WORKERS" --resume 2>&1 | tee -a "$OUT/logs/profiles.log"
"$PYTHON" -u -m aidd_agent.joint_spatial_blocks \
  --profiles "$OUT/profiles" --receipt "$ROOT/properties-validation.json" \
  --output "$OUT/blocks" --resume 2>&1 | tee -a "$OUT/logs/blocks.log"
echo "E096 completed. Review: $OUT/blocks/report.json"
