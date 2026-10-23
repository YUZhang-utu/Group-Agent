#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
MODE=${1:-pilot}
if [[ $# -gt 0 ]]; then shift; fi
OUT=${AIDD_KRAS_OUTPUT:-$HOME/kras-e116}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
case "$MODE" in
  pilot)
    mkdir -p -- "$OUT"
    "${AIDD_PY:-python}" "$ROOT/scripts/run_kras_weekend.py" \
      --output "$OUT/pilot" --sample 1000 --chunk-size 32 --workers 20 --hours 1 "$@" \
      2>&1 | tee -a "$OUT/pilot.log"
    ;;
  full)
    "${AIDD_PY:-python}" - "$OUT/pilot/report.json" <<'PY'
import json, sys
from pathlib import Path
path = Path(sys.argv[1])
if not path.is_file(): raise SystemExit('Run the pilot first and inspect its report.')
report = json.loads(path.read_text())
if report.get('status') != 'complete' or report.get('scope') != 'spread_pilot':
    raise SystemExit('Pilot did not complete; inspect its log before full scanning.')
print('Pilot dual-target full-library linear estimate (hours):', report.get('dual_full_hours_linear_estimate'))
print('This is not a deadline guarantee. Inspect target pass fractions before launch.')
PY
    "${AIDD_PY:-python}" "$ROOT/scripts/run_kras_weekend.py" \
      --output "$OUT/full" --sample 0 --chunk-size 128 --workers 20 --hours 60 \
      --pilot-report "$OUT/pilot/report.json" "$@" \
      2>&1 | tee -a "$OUT/full.log"
    ;;
  *) echo 'Usage: bash scripts/run_kras_weekend.sh pilot|full [runner options]' >&2; exit 2 ;;
esac
