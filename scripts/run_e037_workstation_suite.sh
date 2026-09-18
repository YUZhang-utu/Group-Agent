#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${AIDD_PY:-python}" -u -m aidd_agent.workstation_suite \
  --batch "${AIDD_BATCH:-/mnt/local/hand/yuzhang/aidd/library-precompute-20260911}" \
  --e034 "${E037_E034_SOURCE:-/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118}" \
  --output "${E037_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e037-workstation-suite/run-v1}" \
  --workers "${E037_WORKERS:-16}" --repeats "${E037_REPEATS:-2}" --profile "$@"
