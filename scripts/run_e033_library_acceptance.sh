#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
AIDD_PY=${AIDD_PY:-python}
BATCH=${AIDD_BATCH:-/mnt/local/hand/yuzhang/aidd/library-precompute-20260911}
OUTPUT=${E033_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e033-library-acceptance-$(date +%Y%m%d-%H%M%S)}
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=${E033_THREADS:-20}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

if [[ ! -f "$BATCH/COMPLETE.json" ]]; then
  echo "E032 is incomplete: missing $BATCH/COMPLETE.json." >&2
  echo "Registration is not index completion. Resume the original precompute command in to_human/E032_CONFORMER_CONFLICT.md." >&2
  exit 66
fi

echo "E033 output: $OUTPUT"
echo "Full integrity checking and exact-reference scanning may take time; they are not online query latency."
"$AIDD_PY" -m aidd_agent.library_acceptance \
  --batch "$BATCH" --output "$OUTPUT" --threads "$OMP_NUM_THREADS" "$@"
echo "Return these files: $OUTPUT/report.md and $OUTPUT/report.json (or metrics.csv)."
