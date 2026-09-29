#!/usr/bin/env bash
# Run after the complete E079 audit; all outputs must be fresh.
set -euo pipefail
if [[ $# -ne 2 && $# -ne 4 ]]; then
  echo 'Usage: bash scripts/run_e080_macrocycle_blocks.sh AUDIT_DIR NEW_OUTPUT_ROOT [REGISTRY_SQLITE LIBRARY_ID]' >&2
  exit 2
fi
PYTHON=${AIDD_PYTHON:-python}
AUDIT=$1
OUTPUT=$2
if [[ -e "$OUTPUT" ]]; then
  echo 'Output root must not exist' >&2
  exit 2
fi
registry_args=()
if [[ $# -eq 4 ]]; then
  registry_args=(--registry "$3" --library-id "$4")
fi
# Compare descriptor variants at one capacity with mandatory validation per variant.
for variant in backbone chemistry typed; do
  "$PYTHON" -m aidd_agent.verified_macrocycle_blocks \
    --audit "$AUDIT" --output "$OUTPUT/$variant" --variant "$variant" \
    --capacities 20000 "${registry_args[@]}"
  "$PYTHON" -m aidd_agent.conformer_block_validate \
    --build "$OUTPUT/$variant" --output "$OUTPUT/validation-$variant" --source-samples 200 \
    --allow-review-for-offline
done
