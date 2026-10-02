#!/usr/bin/env bash
set -euo pipefail
# The output is a new parent directory; sampling artifacts live under samples/.
PROFILE="${1:?Supply sampling.json}"
OUT="${2:?Supply a new E097 evaluation directory}"
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$OUT"
exec 9>"$OUT/run.lock"
flock -n 9 || { echo "Another E097 sampler holds this output lock" >&2; exit 1; }
"${AIDD_PYTHON:-python}" -u -m aidd_agent.block_sampling --profile "$PROFILE" \
  --output "$OUT/samples" --schemes E094 E095 E096 --count 100 --seed 20261002 --resume \
  2>&1 | tee -a "$OUT/sampling.log"
