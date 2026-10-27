#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
INSTALL=${AF3_INSTALL:-/mnt/medchem_taltio/wrk/yu_agent/alphafold3}
IMAGE=${AF3_IMAGE:-$INSTALL/images/alphafold3-local-rebuilt.sif}
MODELS=${AF3_MODELS:-$INSTALL/models}
DB_ROOT=/mnt/local/hand/yuzhang/alphafold3
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo 'Usage: bash run_af3.sh INPUT.json [--check]' >&2; exit 2
fi
INPUT=$1
MODE=${2:-run}
if [[ "$MODE" != run && "$MODE" != --check && "$MODE" != --data-only && "$MODE" != --inference-only ]]; then
  echo 'Usage: bash run_af3.sh INPUT.json [--check|--data-only|--inference-only]' >&2; exit 2
fi
[[ -f "$INPUT" && -r "$INPUT" ]] || { echo "Unreadable input: $INPUT" >&2; exit 1; }
INPUT_DIR=$(cd -- "$(dirname -- "$INPUT")" && pwd)
INPUT_NAME=$(basename -- "$INPUT")
INPUT="$INPUT_DIR/$INPUT_NAME"
JOB=${INPUT_NAME%.json}
if [[ "$INPUT_NAME" != *.json || ! "$JOB" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]]; then
  echo 'Use a .json filename with letters, digits, underscores or hyphens.' >&2; exit 2
fi
command -v apptainer >/dev/null || { echo 'Load the school Apptainer module first.' >&2; exit 1; }
[[ -r "$IMAGE" && -d "$MODELS" && -r "$INPUT" ]] || {
  echo "Missing image, model directory or input: $IMAGE | $MODELS | $INPUT" >&2; exit 1;
}
if [[ -n "${AF3_DB_DIR:-}" ]]; then
  DATABASES=$AF3_DB_DIR
else
  candidates=()
  for location in "$DB_ROOT" "$DB_ROOT/databases"; do
    if [[ -d "$location" ]] && compgen -G "$location/uniref90*" >/dev/null; then
      candidates+=("$location")
    fi
  done
  if [[ ${#candidates[@]} -ne 1 ]]; then
    echo 'Database root is absent or ambiguous. Set AF3_DB_DIR to the directory containing the actual AF3 database files.' >&2
    printf 'Checked: %s\n' "$DB_ROOT" "$DB_ROOT/databases" >&2; exit 1
  fi
  DATABASES=${candidates[0]}
fi
[[ -d "$DATABASES" && -r "$DATABASES" ]] || { echo "Unreadable database directory: $DATABASES" >&2; exit 1; }
CACHE=${AF3_JAX_CACHE_DIR:-$SCRIPT_DIR/.af3_cache/jax}
mkdir -p -- "$CACHE"
CACHE=$(cd -- "$CACHE" && pwd)
for location in "$INPUT_DIR" "$IMAGE" "$MODELS" "$DATABASES" "$CACHE"; do
  case "$location" in *:*|*,*|*$'\n'*) echo "Unsupported bind-path separator: $location" >&2; exit 1;; esac
done
printf 'Input: %s\nImage: %s\nModels: %s\nDatabases: %s\n' "$INPUT" "$IMAGE" "$MODELS" "$DATABASES"
printf 'JAX persistent cache: %s\nStage: %s\n' "$CACHE" "$MODE"
echo 'Database file names (top level):'
find "$DATABASES" -mindepth 1 -maxdepth 1 -printf '%f\n' | sort
if [[ "$MODE" == --check ]]; then
  echo 'Host path checks passed. GPU allocation, complete database contents and container compatibility are not yet tested.'
  exit 0
fi
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"
OUTPUT="$INSTALL/tests/output/$JOB/$RUN_ID"
LOG="$INSTALL/logs/${JOB}_$RUN_ID.log"
mkdir -p -- "$OUTPUT" "$(dirname -- "$LOG")"
cp -- "$INPUT" "$OUTPUT/submitted_input.json"
GPU_ARGS=(--nv)
STAGE_ARGS=(--run_data_pipeline=true --run_inference=true)
if [[ "$MODE" == --data-only ]]; then
  GPU_ARGS=()
  STAGE_ARGS=(--run_data_pipeline=true --run_inference=false)
elif [[ "$MODE" == --inference-only ]]; then
  STAGE_ARGS=(--run_data_pipeline=false --run_inference=true)
fi
COMMAND=(apptainer exec "${GPU_ARGS[@]}"
  --bind "$INPUT_DIR:/root/af_input:ro"
  --bind "$OUTPUT:/root/af_output:rw"
  --bind "$MODELS:/root/models:ro"
  --bind "$DATABASES:/root/public_databases:ro"
  --bind "$CACHE:/root/af_jax_cache:rw"
  "$IMAGE" python /app/alphafold/run_alphafold.py
  "--json_path=/root/af_input/$INPUT_NAME"
  --model_dir=/root/models --db_dir=/root/public_databases
  --output_dir=/root/af_output --flash_attention_implementation=xla
  --jax_compilation_cache_dir=/root/af_jax_cache "${STAGE_ARGS[@]}")
{
  printf 'Output: %s\nLog: %s\nCommand: ' "$OUTPUT" "$LOG"
  printf '%q ' "${COMMAND[@]}"
  printf '\n'
} | tee "$LOG"
"${COMMAND[@]}" 2>&1 | tee -a "$LOG"
