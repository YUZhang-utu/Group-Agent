#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export DGLBACKEND=pytorch OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# AIDD_PY is the existing AIDD Python 3.11+ interpreter. The profile independently
# selects the isolated Python 3.9 interpreter for EquiScore GPU inference.
"${AIDD_PY:-python}" -u -m aidd_agent.equiscore "$@"
