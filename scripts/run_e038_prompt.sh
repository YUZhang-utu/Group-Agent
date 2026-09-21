#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ARGS=(--storage-root "${AIDD_PROMPT_STORAGE:-/mnt/local/hand/yuzhang/aidd/prompt-workspace}")
if [[ -n "${AIDD_RUNTIME_PROFILE:-}" ]]; then ARGS+=(--runtime "$AIDD_RUNTIME_PROFILE"); fi
"${AIDD_PY:-python}" "$ROOT/scripts/prompt_aidd.py" "${ARGS[@]}" "$@"
