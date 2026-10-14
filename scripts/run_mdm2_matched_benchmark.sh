#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
OUTPUT="${MDM2_BENCHMARK_OUT:-$HOME/mdm2-e109}"
PYTHON="${PYTHON:-python}"
STAGE="${1:-help}"
WORKERS="${WORKERS:-20}"
FRACTION="${FRACTION:-0.01}"
TOP_BLOCKS="${TOP_BLOCKS:-20}"
REFERENCE_K="${REFERENCE_K:-1000}"
cli=("$PYTHON" scripts/benchmark_matched_blocks.py)
common=(--prepared "$OUTPUT/prepared" --workers "$WORKERS" --chunk-size 4048 --output-budget 100000)
case "$STAGE" in
  prepare)
    extra=()
    if [[ -n "${MDM2_CAMPAIGN:-}" ]]; then extra+=(--campaign "$MDM2_CAMPAIGN"); fi
    "${cli[@]}" prepare --output "$OUTPUT/prepared" "${extra[@]}"
    ;;
  pilot)
    : "${HARDWARE:?Set HARDWARE to your CPU/RAM/storage description}"
    "${cli[@]}" evaluate "${common[@]}" --hardware "$HARDWARE" --fraction "$FRACTION" --mode full --output "$OUTPUT/pilot-full"
    "${cli[@]}" evaluate "${common[@]}" --hardware "$HARDWARE" --fraction "$FRACTION" --mode blocks --scheme E094 --scorer chemplp --top-blocks all --output "$OUTPUT/pilot-all"
    "${cli[@]}" verify-equivalence --full "$OUTPUT/pilot-full" --blocks "$OUTPUT/pilot-all" --output "$OUTPUT/pilot-equivalence"
    "${cli[@]}" replay --prepared "$OUTPUT/prepared" --reference "$OUTPUT/pilot-full" --fractions "$FRACTION" --top-k 100 1000 10000 --rerank-exports --output "$OUTPUT/pilot-replay"
    ;;
  full)
    : "${HARDWARE:?Set HARDWARE to your CPU/RAM/storage description}"
    test -f "$OUTPUT/pilot-equivalence/report.json"
    "${cli[@]}" evaluate "${common[@]}" --hardware "$HARDWARE" --fraction 1 --mode full --resume --output "$OUTPUT/full"
    ;;
  replay)
    "${cli[@]}" replay --prepared "$OUTPUT/prepared" --reference "$OUTPUT/full" --fractions 0.01 0.05 0.1 0.25 0.5 1 --top-k 100 1000 10000 --rerank-exports --output "$OUTPUT/replay"
    ;;
  timing)
    : "${HARDWARE:?Set HARDWARE to your CPU/RAM/storage description}"
    test -f "$OUTPUT/pilot-equivalence/report.json"
    methods=(full E094/chemplp E094/equiscore E095/chemplp E095/equiscore E096/chemplp E096/equiscore)
    runs=()
    for repeat in 1 2 3; do
      for offset in 0 1 2 3 4 5 6; do
        method="${methods[$(((offset + repeat - 1) % 7))]}"
        destination="$OUTPUT/timing-$FRACTION-top$TOP_BLOCKS/$repeat/${method//\//-}"
        mode=(--mode full)
        if [[ "$method" != full ]]; then mode=(--mode blocks --scheme "${method%/*}" --scorer "${method#*/}" --top-blocks "$TOP_BLOCKS"); fi
        "${cli[@]}" evaluate "${common[@]}" --hardware "$HARDWARE" --fraction "$FRACTION" --repeat "$repeat" "${mode[@]}" --output "$destination"
        runs+=("$destination")
      done
    done
    "${cli[@]}" timing-summary --runs "${runs[@]}" --reference-k "$REFERENCE_K" --target-recall 0.95 --output "$OUTPUT/timing-$FRACTION-top$TOP_BLOCKS/summary"
    ;;
  *)
    printf '%s\n' 'Stages: prepare, pilot, full, replay, timing (run separately in that order).' \
      'Set HARDWARE before evaluation; activate the existing AIDD scientific environment.' \
      'Optional: MDM2_BENCHMARK_OUT, MDM2_CAMPAIGN, WORKERS, FRACTION, TOP_BLOCKS, REFERENCE_K.' \
      'Timing runs 3 fresh full and 18 fresh prefix searches; it is not a cheap replay.' \
      'Read to_human/E109_MATCHED_BLOCK_BENCHMARK.md before full evaluation.'
    ;;
esac
