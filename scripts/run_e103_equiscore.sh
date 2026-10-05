#!/usr/bin/env bash
# Workstation defaults for the user's completed, terminal-submitted PLANTS panel.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE=${PLANTS_REPORT:-/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/plants-parallel20/report.json}
PROFILE=${EQUISCORE_PROFILE:-/mnt/local/hand/yuzhang/aidd/tools/equiscore/profile.json}
OUTPUT=${EQUISCORE_OUTPUT:-$(dirname -- "$(dirname -- "$SOURCE")")}
case ${1:---help} in
  pilot)
    bash "$ROOT/scripts/run_equiscore.sh" run --source "$SOURCE" --profile "$PROFILE" --output "$OUTPUT/equiscore-pilot" --scope pilot
    ;;
  full)
    bash "$ROOT/scripts/run_equiscore.sh" run --source "$SOURCE" --profile "$PROFILE" --output "$OUTPUT/equiscore-full" --scope full --pilot-report "$OUTPUT/equiscore-pilot/report.json"
    ;;
  analyze)
    bash "$ROOT/scripts/run_equiscore.sh" analyze --source "$OUTPUT/equiscore-full/report.json" --output "$OUTPUT/equiscore-analysis"
    ;;
  *)
    echo 'Usage: bash scripts/run_e103_equiscore.sh pilot|full|analyze'
    echo 'Run from the existing AIDD environment (Python 3.11+).'
    echo 'Overrides: AIDD_PY, PLANTS_REPORT, EQUISCORE_PROFILE, EQUISCORE_OUTPUT.'
    [[ ${1:---help} == --help ]]
    ;;
esac
