#!/usr/bin/env bash
set -euo pipefail

# Recreate evaluation evidence only. Never rebuild the accepted library or
# repair/overwrite old markers. E034 runs only after E033 exits successfully.
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
STAMP=$(date +%Y%m%d-%H%M%S)-$$
export E033_OUTPUT="${E034_RECOVERY_ROOT:-/mnt/local/hand/yuzhang/aidd}/e033-library-acceptance/recheck-$STAMP"
export E034_OUTPUT="${E034_RECOVERY_ROOT:-/mnt/local/hand/yuzhang/aidd}/e034-expanded-wee1/recheck-$STAMP"
echo "Fresh E033 evidence: $E033_OUTPUT"
echo "Fresh E034 results: $E034_OUTPUT"
bash "$ROOT/scripts/run_e033_library_acceptance.sh"
bash "$ROOT/scripts/run_e034_expanded_wee1.sh"
