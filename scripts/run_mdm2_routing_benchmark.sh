#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
REFERENCE="${MDM2_REFERENCE:-/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers}"
RUNS="${MDM2_RUNS_ROOT:-/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs}"
OUTPUT="${1:-$HOME/mdm2-routing-benchmark-$(date +%Y%m%d-%H%M%S)}"
if [[ $# -gt 0 ]]; then shift; fi
python scripts/benchmark_routing.py historical \
  --reference "$REFERENCE" --runs-root "$RUNS" --output "$OUTPUT" "$@"
