# E038: prompts for protein data, AF3 and calibrated 3D search

For the user-supplied Apptainer installation and selectable GPT/DeepSeek profiles,
use [the container/provider setup](AF3_CONTAINER_AND_LLM_PROVIDERS.md) instead of
the native AF3 configuration below.

The interface uses OpenAI-compatible Chat Completions. The user reports AF3 already
installed on the workstation. The LLM generates a structured plan; local adapters
retrieve verified protein sequences, prepare AF3 JSON and invoke the configured
runner. AF3 consumes structured input, not a free-text prompt.

All maintained instructions, examples and planner summaries are English. Repository
skills are in `.agents/skills/`; see [the skill guide](../docs/workflow-skills.md).

## Workstation sequence

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main

# No network, API key or GPU required.
bash scripts/run_e038_prompt_smoke.sh

# Real-library performance/equivalence suite, if still pending.
bash scripts/run_e037_workstation_suite.sh
```

The smoke report under `e038-prompt-smoke/<timestamp>-<PID>/report.json` is explicitly
synthetic and does not validate live LLM/protein APIs or AF3. E037 measures the real
library under its separate protocol.

## Configure live services

Store local configuration outside the repository, for example under
`/mnt/local/hand/yuzhang/aidd/config/`. Copy
`configs/llm/openai-compatible.example.json`; set the actual `base_url` (through
`/v1`) and `model`. Keep `api_key_env` as `AIDD_LLM_API_KEY`. Set the secret in the
terminal environment, not JSON, Git or conversation messages.

A local compatible server may use `http://127.0.0.1:<port>/v1`; omit `api_key_env`
when authentication is unnecessary. If JSON mode is unsupported, set `json_mode`
to false; strict local response validation still applies.

```bash
export AIDD_LLM_PROFILE=/mnt/local/hand/yuzhang/aidd/config/llm.local.json
read -rsp 'LLM API Key: ' AIDD_LLM_API_KEY
echo
export AIDD_LLM_API_KEY
```

Copy `configs/prompt-runtime.example.json` outside the repository and point its
`af3_profile` at the existing local AF3 profile. See
`configs/models/alphafold3.example.json`: `python`, `runner`, `model_parameters` and
`databases` must be real absolute paths and satisfy the existing AF3 profile checks.
This tool does not install AF3, download weights or alter its environment. The AF3
pipeline handles MSA/templates; the LLM does not fabricate them. Record the actual
installed `version` in the profile when available.

```bash
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime.local.json
```

The example search profile retains reference seed generation
(`bounded_pair_seeds:false`). After E037 passes, choose scheduling/engine settings
from its results. Changed configuration requires a new plan, not mixed old receipts.

## Test the live planner, then prepare inputs

```bash
# LLM only: no protein API calls or scientific execution.
bash scripts/run_e038_prompt.sh --plan-only \
  --prompt 'Find human WEE1 in UniProt and PDB candidates at resolution at most 3 Angstrom. Prepare full-length AF3 input but do not run prediction.'

# Retrieve evidence and prepare FASTA, provenance, PDB candidates and AF3 JSON.
bash scripts/run_e038_prompt.sh \
  --prompt 'Find human WEE1 in UniProt and PDB candidates at resolution at most 3 Angstrom. Prepare full-length AF3 input but do not run prediction.'
```

The wrapper creates/reuses and activates the `workstation` user's `Prompt AIDD`
Project. Default storage is `/mnt/local/hand/yuzhang/aidd/prompt-workspace`;
override with `AIDD_PROMPT_STORAGE`, `--username` or `--project-name`. Existing
Projects can use the low-level CLI with explicit `--db --user --project`; ownership
and activation are checked.

Printed output includes `plan.json`, context and report paths. Tasks live under
Project `runs/PROMPT-<id>/`; requests/plans are also recorded in the existing AI audit
tables. Plans allow fixed actions, not shell commands, arbitrary paths, invented
sequences or model-chosen budgets. Cloud context contains only the user prompt and
capability/skill descriptions; local library, sequence and experimental data are
not automatically uploaded.

## Run installed AF3

```bash
bash scripts/run_e038_prompt.sh --allow-compute \
  --prompt 'Retrieve the full-length human WEE1 UniProt sequence, prepare AF3 input and run one AF3 prediction with seed 1.'
```

You can explicitly request an accession, a1-based inclusive residue range or aCCD
ligand: 'Retrieve full-length P30291 and prepare an AF3 complex with ATP; do not run
prediction yet.' CCD identity is checked against RCSB. The current prompt adapter
supports one protein chain plus CCD ligands, not DNA/RNA, multiple proteins, custom
CCD/covalent bonds or automatic domain selection. Missing information produces
clarification questions. Multiple UniProt matches are not silently resolved, and
nonstandard residues stop preparation rather than being replaced.

`--allow-compute` enables AF3/3D execution for that invocation. Without it, data
retrieval/preparation can complete and a compute step reports `blocked`. Model
instructions cannot bypass this flag. Successful AF3 execution records a model,
confidence and `prediction-manifest.json` with hashes; it does not establish pose
or biological quality. Failed attempts are preserved in separate directories.

## Prompt-driven 3D search

```bash
bash scripts/run_e038_prompt.sh --allow-compute \
  --prompt 'Run the calibrated WEE1 QT9 query against the existing full library with Gaussian refinement. Keep E031 annotation-only.'
```

Available queries: `wee1_qt9`, `wee1_824`, `wee1_both`. Another target is never silently
mapped to WEE1. New targets still require prepared and validated queries/receptors.
The model cannot change the10000candidate budget, Top-N or seed cap, nor promote
E031 to filtering/ranking. Search preserves the original E034 equivalence checks.

## Resume and return evidence

Use the printed real IDs/paths to resume the same plan:

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m aidd_agent.prompt_workflow run \
  --db /actual/registry/aidd.sqlite3 --user USR-ACTUAL --project PRJ-ACTUAL \
  --plan /actual/Project/runs/PROMPT-ACTUAL/plan.json \
  --runtime "$AIDD_RUNTIME_PROFILE" --allow-compute
```

Plan/code/configuration/output hashes must match. `ask` and the wrapper create new
plans; they do not discover old plans automatically. Reports distinguish complete,
blocked, failed and not_run; inspect the failed step's local execution log.
Return the smoke report, E037 summary, and live task `plan.json`/
`execution/report.json`, excluding credentials and sensitive prompt text.

## References and validation limits

- [OpenAI Chat API](https://developers.openai.com/api/reference/resources/chat)
- [UniProt REST queries](https://www.uniprot.org/help/api_queries)
- [RCSB Data API](https://data.rcsb.org/)
- [Official AF3 input format](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)

Basic AF3 dialect-v1 protein/CCD inputs support existing installations. Historical
local validation:222tests passed/2dependency skips and a synthetic smoke run.
Live LLM, protein APIs and AF3/GPU execution were not tested locally; mock checks
do not replace workstation integration validation.
