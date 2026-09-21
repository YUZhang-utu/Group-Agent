#!/usr/bin/env bash
set -euo pipefail
# Run after the operator's normal module load. This never runs a docking job.
if [[ $# -gt 1 ]]; then printf 'Usage: bash scripts/probe_docking_tools.sh [module-name]\n' >&2; exit 2; fi
if [[ $# -eq 1 ]]; then
  if ! type module >/dev/null 2>&1; then printf 'The module command is unavailable in this shell. Load your site environment first.\n' >&2; exit 2; fi
  module load "$1"
fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
PROFILE="${AIDD_DOCKING_PROFILE:-${AIDD_CONFIG_DIR:-/mnt/local/hand/yuzhang/aidd/config}/docking.local.json}"
ARGS=()
if [[ -n "${AIDD_DOCKING_PROFILE:-}" || -f "$PROFILE" ]]; then ARGS+=(--profile "$PROFILE"); fi
"${AIDD_PY:-python}" -m aidd_agent.docking_discovery "${ARGS[@]}"
