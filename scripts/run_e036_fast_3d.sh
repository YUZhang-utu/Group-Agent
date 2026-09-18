#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${AIDD_PY:-python}" -m aidd_agent.fast_3d_search \
  --batch "${AIDD_BATCH:-/mnt/local/hand/yuzhang/aidd/library-precompute-20260911}" \
  --e034 "${E036_E034_SOURCE:-/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118}" \
  --output "${E036_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e036-fast-3d/run-w16-c500-r64-v1}" \
  --workers "${E036_WORKERS:-16}" --coarse-chunk "${E036_COARSE_CHUNK:-500}" \
  --refine-chunk "${E036_REFINE_CHUNK:-64}" --fresh-retrieval "$@"
