#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG_DIR="${AIDD_CONFIG_DIR:-/mnt/local/hand/yuzhang/aidd/config}"
mkdir -p "$CONFIG_DIR/llm"
copy_new() {
  if [[ -e "$2" ]]; then
    printf 'Preserved existing profile: %s\n' "$2"
  else
    cp -- "$1" "$2"
    printf 'Created profile: %s\n' "$2"
  fi
}
copy_new "$ROOT/configs/llm/gpt.json" "$CONFIG_DIR/llm/gpt.json"
copy_new "$ROOT/configs/llm/deepseek.json" "$CONFIG_DIR/llm/deepseek.json"
copy_new "$ROOT/configs/models/alphafold3-apptainer.example.json" "$CONFIG_DIR/alphafold3-apptainer.local.json"
copy_new "$ROOT/configs/prompt-runtime-apptainer.example.json" "$CONFIG_DIR/prompt-runtime-apptainer.local.json"
copy_new "$ROOT/configs/docking/workstation.example.json" "$CONFIG_DIR/docking.local.json"
printf 'Set AIDD_LLM_CONFIG_DIR=%s/llm\n' "$CONFIG_DIR"
printf 'Set AIDD_RUNTIME_PROFILE=%s/prompt-runtime-apptainer.local.json\n' "$CONFIG_DIR"
printf 'Confirm AF3 terms and set license_acknowledged=true in your local AF3 profile before prediction.\n'
