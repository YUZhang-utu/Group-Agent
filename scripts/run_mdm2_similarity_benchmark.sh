#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
SOURCE="${MDM2_BENCHMARK_OUT:-$HOME/mdm2-e109}"
OUTPUT="${SIMILARITY_BENCHMARK_OUT:-$HOME/mdm2-e110-e094}"
PYTHON="${PYTHON:-python}"
STAGE="${1:-help}"
WORKERS="${WORKERS:-20}"
SCHEMES="${SCHEMES:-E094}"
read -r -a schemes <<< "$SCHEMES"
cli=("$PYTHON" scripts/benchmark_similarity_blocks.py)
plan() {
  "${cli[@]}" plan --prepared "$SOURCE/prepared" --reference "$SOURCE/pilot-full" --schemes "${schemes[@]}" --sample-size 500 --output "$OUTPUT/plan"
}
fit() {
  : "${HARDWARE:?Set HARDWARE to the CPU/RAM/storage description}"
  "${cli[@]}" fit --plan "$OUTPUT/plan" --workers "$WORKERS" --hardware "$HARDWARE" --resume --output "$OUTPUT/fit"
}
replay() {
  "${cli[@]}" replay --plan "$OUTPUT/plan" --fit "$OUTPUT/fit" --output "$OUTPUT/replay"
}
bundle() {
  tar -czf "$OUTPUT/reports.tar.gz" -C "$OUTPUT" plan/report.json plan/protocol.json fit/report.json fit/priorities.csv replay
  printf "Reports: %s\n" "$OUTPUT/reports.tar.gz"
}
case "$STAGE" in
  plan) plan ;;
  fit) fit ;;
  replay) replay; bundle ;;
  all) : "${HARDWARE:?Set HARDWARE to the CPU/RAM/storage description}"; plan; fit; replay; bundle ;;
  bundle) bundle ;;
  *) printf "%s\n" "Stages: plan, fit, replay, all, bundle" "Default: E094, exactly 500 distinct molecules per block, existing E109 pilot validation." "Plan prints missing work and estimated hours; fit computes only missing sample conformers." "Use a fresh output for all; after an interruption resume with fit, then replay." "Optional: SCHEMES=E094 E095 E096, SIMILARITY_BENCHMARK_OUT, MDM2_BENCHMARK_OUT, WORKERS, PYTHON." ;;
esac
