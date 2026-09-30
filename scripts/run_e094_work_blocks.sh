#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/run_e094_work_blocks.sh base|properties|all BACKBONE_BUILD OUTPUT_ROOT
# Keep OUTPUT_ROOT outside the source/model/routing directories. No old inputs are rewritten.
MODE="${1:-base}"
BUILD="${2:?Supply the completed backbone build directory}"
OUT="${3:?Supply a new E094 output root}"
PYTHON="${AIDD_PYTHON:-python}"
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT/logs" "$OUT/sqlite-tmp"
export SQLITE_TMPDIR="$OUT/sqlite-tmp"

if [[ "$MODE" != base && "$MODE" != properties && "$MODE" != all ]]; then
  echo "Mode must be base, properties or all" >&2
  exit 2
fi

if [[ "$MODE" == base || "$MODE" == all ]]; then
  "$PYTHON" -u -m aidd_agent.final_work_blocks \
    --model "$BUILD/capacity-20000" \
    --diagnostic "${BUILD}-fragmentation-review" \
    --candidates "${BUILD}-boundary-candidates-v1" \
    --routing "${BUILD}-calibrated-routing-v1" \
    --output "$OUT/backbone" --minimum-block-size "${AIDD_MIN_BLOCK_SIZE:-5000}" \
    2>&1 | tee "$OUT/logs/backbone.log"
  "$PYTHON" -u -m aidd_agent.work_block_delivery validate \
    --blocks "$OUT/backbone" --output "$OUT/backbone-validation.json" \
    2>&1 | tee "$OUT/logs/backbone-validation.log"
fi

if [[ "$MODE" == properties || "$MODE" == all ]]; then
  "$PYTHON" -u -m aidd_agent.work_block_properties prepare \
    --blocks "$OUT/backbone" --build "$BUILD" --output "$OUT/profiles" \
    --workers "${AIDD_PROPERTY_WORKERS:-8}" --resume \
    2>&1 | tee -a "$OUT/logs/profiles.log"
  "$PYTHON" -u -m aidd_agent.work_block_properties refine \
    --profiles "$OUT/profiles" --output "$OUT/properties" \
    --maximum-children 8 --minimum-contrast 0.35 \
    2>&1 | tee "$OUT/logs/properties.log"
  "$PYTHON" -u -m aidd_agent.work_block_delivery validate \
    --blocks "$OUT/backbone" --properties "$OUT/properties" \
    --output "$OUT/properties-validation.json" \
    2>&1 | tee "$OUT/logs/properties-validation.log"
fi

echo "E094 requested stages complete. Review reports and validation receipts under: $OUT"
