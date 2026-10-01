#!/usr/bin/env bash
set -euo pipefail
# Run from the Group-Agent checkout; E094 sources stay unchanged.
ROOT="${1:?Supply the completed final-blocks-e094 directory}"
OUT="${2:?Supply a separate E095 output directory}"
PYTHON="${AIDD_PYTHON:-python}"
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
"$PYTHON" -u -m aidd_agent.bounded_property_blocks \
  --profiles "$ROOT/profiles" --receipt "$ROOT/properties-validation.json" \
  --output "$OUT" --resume 2>&1 | tee -a "${OUT}.log"
echo "E095 complete. Review: $OUT/report.json and $OUT/split-decisions.json"
