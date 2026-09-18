#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
AIDD_PY=${AIDD_PY:-python}
BATCH=${AIDD_BATCH:-/mnt/local/hand/yuzhang/aidd/library-precompute-20260911}
E033=${E033_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604}
OUTPUT=${E034_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/20260918-run1}
QT9=${E034_QT9_DIR:-$ROOT/data/e019_query_8bju}
X8B=${E034_1X8B_DIR:-$ROOT/data/e026_query_1x8b}
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
# Prevent BLAS oversubscription in Gaussian worker processes; FAISS threads
# are independently configured by --threads in the Python driver.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "E033 evidence: $E033"
echo "E034 output: $OUTPUT"
echo "Preflight checks and exact reference scans are not online query latency."
"$AIDD_PY" -m aidd_agent.expanded_wee1 \
  --batch "$BATCH" --e033 "$E033" --output "$OUTPUT" \
  --qt9 "$QT9" --x8b "$X8B" \
  --workers "${E034_WORKERS:-16}" --threads "${E034_THREADS:-20}" "$@"
echo "Return: $OUTPUT/report.md and $OUTPUT/report.json"
