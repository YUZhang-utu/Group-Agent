# Installed AF3 container and selectable GPT / DeepSeek profiles

This setup preserves the supplied Apptainer installation. No container rebuild,
AF3 installation or model download is performed. Native Python AF3 profiles remain
supported. All authored configuration and guidance is English.

## Set up local profiles

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/setup_prompt_profiles.sh
export AIDD_LLM_CONFIG_DIR=/mnt/local/hand/yuzhang/aidd/config/llm
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime-apptainer.local.json
```

The setup script preserves existing files. It copies these editable local profiles:

- `llm/gpt.json`: `https://api.openai.com/v1`, `gpt-5.6-terra`, `OPENAI_API_KEY`.
- `llm/deepseek.json`: `https://api.deepseek.com`, `deepseek-flash`, `DEEPSEEK_API_KEY`.
- `alphafold3-apptainer.local.json`: your supplied SIF, weights and database paths.
- `prompt-runtime-apptainer.local.json`: the relative AF3 profile and existing search inputs.

Set `AIDD_CONFIG_DIR` before setup to choose another destination. Set the exported
LLM/runtime paths to match that destination. Change model IDs in the local JSON
when needed; account availability still needs a real plan-only test.

Only key ENVIRONMENT VARIABLE NAMES belong in JSON. Enter either or both secrets:

```bash
read -rsp 'OpenAI API Key: ' OPENAI_API_KEY
printf '\n'
export OPENAI_API_KEY
read -rsp 'DeepSeek API Key: ' DEEPSEEK_API_KEY
printf '\n'
export DEEPSEEK_API_KEY
```

You only need the selected provider's key. Keys are not shared across providers and
there is no automatic fallback. An explicit `--provider` overrides the legacy
`AIDD_LLM_PROFILE` environment variable. An explicit `--llm-profile` selects a custom
file and cannot be combined with `--provider`. Without either flag, the legacy
profile is honored; otherwise `AIDD_LLM_PROVIDER` selects a provider (default GPT).
`--llm-config-dir` overrides `AIDD_LLM_CONFIG_DIR`; without either, bundled profiles
are used. Missing local profiles are errors, not silent fallback to another service.

## Test each real planner

```bash
bash scripts/run_e038_prompt.sh --provider deepseek --plan-only \
  --prompt 'Retrieve human WEE1 P30291 and prepare full-length AF3 input with seed 1. Do not run prediction.'
bash scripts/run_e038_prompt.sh --provider gpt --plan-only \
  --prompt 'Retrieve human WEE1 P30291 and prepare full-length AF3 input with seed 1. Do not run prediction.'
```

Inspect the printed plans for correct target/actions/dependencies. A valid JSON
response is not by itself evidence of semantic correctness. Switching providers
creates a new plan; resuming a sealed plan does not call another LLM.

## Check the installed AF3 profile

The supplied profile uses:

```text
execution: apptainer
executable: apptainer (resolved from the workstation PATH)
image: /mnt/medchem_taltio/wrk/yu_agent/alphafold3/images/alphafold3-local-rebuilt.sif
python inside container: python
runner inside container: /app/alphafold/run_alphafold.py
models: /mnt/medchem_taltio/wrk/yu_agent/alphafold3/models
databases: /mnt/local/hand/yuzhang/alphafold3/databases
flash_attention_implementation: xla
```

The existing AF3 model-terms check is preserved. After confirming the terms for
this installation, set `license_acknowledged` to `true` in your LOCAL AF3 profile.
The distributed example does not assert that acknowledgement on your behalf.

```bash
command -v apptainer
nvidia-smi
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python - <<'PY'
import json, os
from pathlib import Path
from aidd_agent.prediction import load_model_profile, validate_af3_installation
runtime = Path(os.environ['AIDD_RUNTIME_PROFILE'])
cfg = json.loads(runtime.read_text())
p = Path(cfg['af3_profile'])
if not p.is_absolute():
    p = runtime.parent / p
profile = load_model_profile(p)
validate_af3_installation(profile)
print('Host AF3 prerequisites passed; actual GPU inference still needs validation.')
PY
```

This checks the host SIF/directories/executable, not the container's Python packages
or GPU compatibility. Load your normal Apptainer module first if necessary.

## Execute AF3 using the selected provider

```bash
bash scripts/run_e038_prompt.sh --provider deepseek --allow-compute \
  --prompt 'Retrieve full-length human WEE1 P30291, prepare AF3 input with seed 1 and run one prediction.'
```

Replace `deepseek` with `gpt` to select GPT. Each request uses its own input and
`prediction-attempt-NNN` output directory. Inputs, weights and databases are mounted
read-only; the output mount is writable. The command retains `exec --nv` and XLA,
uses a list of arguments without a shell pipeline, and captures stdout/stderr in
`execution.log`. Nonzero container exit codes fail the step. A valid final model
and confidence output are required before reporting `prediction_completed`.
`command.json` records the actual argv. Failed outputs are retained; retries use a
new attempt directory. Completed valid stages are reused through the existing
low-level resume command in the full acceptance guide.

Image SHA-256 is checked once per workflow invocation and recorded in the execution
protocol; changed images invalidate resume. Reading a large SIF is startup integrity
cost and can take time. It is not online 3D-search latency. Model/database trees are
not recursively hashed by this adapter. Preserve their versions separately.

The historical `WEE1_8BJU_puuttuva_nomsa.json` is not reproduced merely from its
filename. Current preparation uses the verified requested sequence; existing
MSA/template/construct settings must not be inferred from that name. Supporting
this container does not add arbitrary AF3 JSON ingestion to the prompt interface.

## Full combined acceptance

Use [the full acceptance guide](WORKSTATION_FULL_ACCEPTANCE.md), replacing its old
native AF3 setup section with this guide. E037 can run before either API is configured.
For the combined prompt in step 5, add `--provider deepseek` or `--provider gpt`.
Do not run AF3 concurrently with E037 performance measurement. Use fresh plans after
code/runtime changes; preserve previous E033-E035 evidence.

Local checks simulate HTTP and container execution. Actual API billing/access,
container runner compatibility and GPU inference require the workstation.

## Official references (checked 2026-09-21)

- [DeepSeek current API setup and model names](https://api-docs.deepseek.com/)
- [DeepSeek JSON output](https://api-docs.deepseek.com/guides/json_mode/)
- [GPT-5.6 Terra](https://developers.openai.com/api/docs/models/gpt-5.6-terra)
- [Apptainer bind mounts](https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html)
