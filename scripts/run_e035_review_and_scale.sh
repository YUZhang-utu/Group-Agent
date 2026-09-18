#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
AIDD_PY=${AIDD_PY:-python}
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
BATCH=${AIDD_BATCH:-/mnt/local/hand/yuzhang/aidd/library-precompute-20260911}
SOURCE=${E035_E034_SOURCE:-/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118}
OUTPUT=${E035_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/run-100k-v1}
echo "E034 input: $SOURCE"
echo "E035 output: $OUTPUT"
echo "Scale: ${E035_SCALE:-100000} distinct library conformers; constructed stress poses, not docking."
"$AIDD_PY" -m aidd_agent.e035_validation --batch "$BATCH" --e034 "$SOURCE" --output "$OUTPUT" \
  --scale "${E035_SCALE:-100000}" --repeats "${E035_REPEATS:-3}" "$@"
echo "Return $OUTPUT/report.md and report.json; pose-review files are under each query/review/."
