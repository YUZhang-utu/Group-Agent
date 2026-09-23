#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export AIDD_ASSIGNMENT_BACKEND=numba
batch=/mnt/local/hand/yuzhang/aidd/library-precompute-20260911
proposal=/mnt/local/hand/yuzhang/aidd/protein-only-check-20260923-132230/proposal/report.json
result_root=${1:?Pass a persistent NEW result directory, or the same directory to resume}
budget=${AIDD_RETRIEVAL_MOLECULES:-1000000}
python -c 'import faiss, numba, rdkit; print("FAISS, Numba and RDKit available")'
mkdir -p "$result_root"
python -m aidd_agent.budget_screen --batch "$batch" --recommendation "$proposal" \
  --definitions to_human/E058_MDM2_DUAL_REGIONS.json --output "$result_root" \
  --retrieval-molecules "$budget" --workers 24 --chunk-conformers 64 \
  --template-quota 100000 --nprobe 128
python -m aidd_agent.budget_export --batch "$batch" --run "$result_root" \
  --output "$result_root/page-000001-100000" --start-rank 1 --count 100000
echo "Ranking: $result_root/ranking.sqlite"
echo "MOL2 and names: $result_root/page-000001-100000"
