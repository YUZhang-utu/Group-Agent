#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
baseline=/mnt/local/hand/yuzhang/aidd/pilot-256-20260923-135254
batch=/mnt/local/hand/yuzhang/aidd/library-precompute-20260911
result_root="/mnt/local/hand/yuzhang/aidd/e058-review-$(date +%Y%m%d-%H%M%S)"
mkdir "$result_root"
echo "Results: $result_root"
python -m aidd_agent.pilot_review --pilot "$baseline" --batch "$batch" \
  --output "$result_root/review" --recover 5
python -m aidd_agent.spatial_consistency --pilot "$baseline" --batch "$batch" \
  --definitions to_human/E058_MDM2_DUAL_REGIONS.json --output "$result_root/spatial"
python -c 'import numba; print("Numba:",numba.__version__)'
export AIDD_ASSIGNMENT_BACKEND=numba
python -m aidd_agent.consensus_pilot \
  --recommendation /mnt/local/hand/yuzhang/aidd/protein-only-check-20260923-132230/proposal/report.json \
  --batch "$batch" --output "$result_root/compiled-pilot" \
  --molecules 256 --workers 24 --chunk-molecules 4 --seconds 600 --seed 20260923
python -m aidd_agent.pilot_review --pilot "$baseline" --batch "$batch" \
  --output "$result_root/comparison" --recover 0 --compare "$result_root/compiled-pilot"
echo "Review: $result_root/review/report.json"
echo "Spatial: $result_root/spatial/report.json"
echo "Comparison: $result_root/comparison/report.json"
