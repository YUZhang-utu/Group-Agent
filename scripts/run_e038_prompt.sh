#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${AIDD_LLM_PROFILE:?Set AIDD_LLM_PROFILE to your local provider JSON (no API key in JSON)}"
ARGS=(--storage-root "${AIDD_PROMPT_STORAGE:-/mnt/local/hand/yuzhang/aidd/prompt-workspace}" \
      --llm-profile "$AIDD_LLM_PROFILE")
if [[ -n "${AIDD_RUNTIME_PROFILE:-}" ]]; then ARGS+=(--runtime "$AIDD_RUNTIME_PROFILE"); fi
"${AIDD_PY:-python}" "$ROOT/scripts/prompt_aidd.py" "${ARGS[@]}" "$@"
