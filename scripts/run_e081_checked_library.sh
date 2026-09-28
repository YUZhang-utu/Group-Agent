#!/usr/bin/env bash
# Full source audit, all descriptor variants at capacity 20000, and validation gates.
set -euo pipefail
if [[ $# -ne 1 && $# -ne 3 ]]; then
  echo 'Usage: bash scripts/run_e081_checked_library.sh NEW_RUN_ROOT [REGISTRY_SQLITE LIBRARY_ID]' >&2
  exit 2
fi
RUN=$1
if [[ -e "$RUN" ]]; then
  echo 'Run root must not exist; interrupted runs need a new root' >&2
  exit 2
fi
bash scripts/run_e079_macrocycle_prepare.sh "$RUN/audit-stage"
if [[ $# -eq 3 ]]; then
  bash scripts/run_e080_macrocycle_blocks.sh "$RUN/audit-stage/audit" "$RUN/blocks-stage" "$2" "$3"
else
  bash scripts/run_e080_macrocycle_blocks.sh "$RUN/audit-stage/audit" "$RUN/blocks-stage"
fi
echo "Offline block validation passed. Inspect reports under $RUN/blocks-stage/validation-* before use."
