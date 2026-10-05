#!/usr/bin/env bash
# Linux x86_64 only. Installs into a dedicated prefix, never the active AIDD env.
set -euo pipefail
PREFIX=${1:-/mnt/local/hand/yuzhang/aidd/tools/equiscore}
if [[ ${1:-} == --help ]]; then
  echo 'Usage: bash scripts/install_equiscore.sh [absolute-install-prefix]'
  echo 'Requires Linux x86_64, Conda, Git, network access and an NVIDIA GPU.'
  exit 0
fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
[[ $# -le 1 && "$PREFIX" == /* && "$PREFIX" != / ]] || { echo 'Provide one absolute dedicated install prefix.' >&2; exit 2; }
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'This environment recipe targets Linux x86_64.' >&2; exit 2; }
command -v conda >/dev/null || { echo 'Conda is not on PATH. Initialize Conda first.' >&2; exit 2; }
command -v git >/dev/null
command -v nvidia-smi >/dev/null || { echo 'NVIDIA driver tools are missing.' >&2; exit 2; }
nvidia-smi --query-gpu=name,driver_version --format=csv
if [[ -d "$PREFIX" && ! -f "$PREFIX/.aidd-equiscore-install" ]]; then
  echo 'Existing unmanaged prefix: choose a new dedicated directory.' >&2; exit 2
fi
mkdir -p -- "$PREFIX"
PREFIX=$(cd -- "$PREFIX" && pwd)
touch "$PREFIX/.aidd-equiscore-install"
COMMIT=8b2a9289cf7d181fa49de6ac6712260e8c500c4a
REPO="$PREFIX/EquiScore"
if [[ ! -d "$REPO" ]]; then
  git clone https://github.com/Intelligent-Drug-Discovery-Lab/EquiScore.git "$REPO"
  git -C "$REPO" checkout --detach "$COMMIT"
fi
[[ $(git -C "$REPO" rev-parse HEAD) == "$COMMIT" ]] || { echo 'Existing EquiScore revision differs; use a new prefix.' >&2; exit 2; }
git -C "$REPO" diff --exit-code HEAD -- '*.py'
if [[ ! -x "$PREFIX/env/bin/python" ]]; then
  conda create -y --prefix "$PREFIX/env" -c conda-forge python=3.9 pip cudatoolkit=11.3
fi
EPY="$PREFIX/env/bin/python"
"$EPY" -m pip --isolated install 'pip==24.3.1' 'setuptools==69.5.1' 'wheel==0.43.0'
"$EPY" -m pip --isolated install 'torch==1.11.0+cu113' --extra-index-url https://download.pytorch.org/whl/cu113
"$EPY" -m pip --isolated install -r "$ROOT/configs/equiscore-requirements.txt" --find-links https://data.dgl.ai/wheels/repo.html
"$EPY" -m pip check
"$EPY" - "$PREFIX" <<'PY'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
weights = root / 'EquiScore/workdir/official_weight/save_model_screen.pt'
assert hashlib.sha256(weights.read_bytes()).hexdigest() == 'd4367bb73686b2363e238abb778fab55e2924458ec1ced561072bd82f711695d', 'Checkpoint hash mismatch'
profile = root / 'profile.json'
expected = dict(python=str(root / 'env/bin/python'), repository=str(root / 'EquiScore'))
if profile.exists():
    current = json.loads(profile.read_text())
    assert all(current.get(k) == v for k, v in expected.items()), 'Existing profile differs; no overwrite'
else:
    profile.write_text(json.dumps(expected, indent=2) + '\n')
PY
"$EPY" -m pip freeze > "$PREFIX/installed-packages.txt"
export DGLBACKEND=pytorch OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="$PREFIX/env/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$EPY" "$ROOT/scripts/equiscore_worker.py" --doctor "$PREFIX/profile.json" | tee "$PREFIX/doctor.log"
echo "Installation checks passed. Profile: $PREFIX/profile.json"
echo 'Next: run the technical pilot on saved PLANTS poses before full scoring.'
