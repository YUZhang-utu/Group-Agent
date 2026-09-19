#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
"${AIDD_PY:-python}" -m aidd_agent.prompt_smoke --output \
  "${E038_OUTPUT:-/mnt/local/hand/yuzhang/aidd/e038-prompt-smoke/$(date +%Y%m%d-%H%M%S)-$$}"
