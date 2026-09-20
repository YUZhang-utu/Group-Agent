# Complete workstation acceptance: search, live prompts, proteins and AF3

## What already passed

The supplied E033 report passed library acceptance/calibration. E034 completed both
WEE1 queries; its original E031 latency gate passed only for 824. E035 subsequently
passed equivalent optimized scoring on one million distinct conformers and six
million comparisons; both historical latency ratios passed. These are user-supplied
workstation results, not a fresh remote artifact audit. Do not repeat these large
runs just to validate skills. E035 exported poses still require human inspection.

E039 local regression passed 224 tests with two dependency skips. Live model/API/AF3
execution and real E037 optimization timings remain pending.

## 1. Update and check the workstation environment

Run from the existing AIDD environment. If it is not active, activate the existing
`envs/aidd-workstation/bin/activate` environment first.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export AIDD_PY="$(command -v python)"
python --version
python -c 'import numpy, rdkit, faiss; print("Search dependencies available")'
python scripts/check_english.py
python -m pytest -q
bash scripts/run_e038_prompt_smoke.sh
```

Inspect skips and failures instead of assuming every workstation matches the local
224-test result. The smoke report must pass, but is synthetic and uses no real
model, protein service or AF3 GPU. Keep its printed output path.

## 2. Validate real-library speed and equivalence

Run this separately from AF3 or other heavy jobs. Use a fresh output, retaining all
previous evidence:

```bash
export E037_OUTPUT="/mnt/local/hand/yuzhang/aidd/e037-workstation-suite/acceptance-$(date +%Y%m%d-%H%M%S)"
printf '%s\n' "$E037_OUTPUT"
bash scripts/run_e037_workstation_suite.sh
cat "$E037_OUTPUT/report.md"
```

The suite performs four scenarios, two rounds and two queries: 16 timed query
executions, followed by a separate two-query profile run. It already covers E036;
no separate E036 run or repeated E035 million-conformer run is required.

Acceptance: root status complete, child candidate/Gaussian/E031 equivalence passed,
and timing rows eligible. Functional completion does not prove acceleration:
inspect `speedup_vs_old_reference` for each query and scenario. A value at or below
1 means no measured speedup. Startup/integrity/index-load costs remain separately
reported; partial chunk reuse cannot establish full-run speedup.

After interruption, retain the SAME output value and run:

```bash
bash scripts/run_e037_workstation_suite.sh --resume
```

Do not regenerate the timestamp when resuming. New code/configuration needs a fresh
output. A partially reused run may finish correctly but need a fresh timing run.

## 3. Configure the real model and installed AF3

```bash
mkdir -p /mnt/local/hand/yuzhang/aidd/config
export AIDD_LLM_PROFILE=/mnt/local/hand/yuzhang/aidd/config/llm.local.json
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime.local.json
export AIDD_PROMPT_STORAGE=/mnt/local/hand/yuzhang/aidd/prompt-workspace
cp -n configs/llm/openai-compatible.example.json "$AIDD_LLM_PROFILE"
cp -n configs/prompt-runtime.example.json "$AIDD_RUNTIME_PROFILE"
cp -n configs/models/alphafold3.example.json /mnt/local/hand/yuzhang/aidd/config/alphafold3.local.json
```

Edit the local files. Existing files are preserved by `cp -n`.

- LLM: actual `base_url` ending in `/v1`, available `model`, and `api_key_env`.
  Official OpenAI uses the example base URL. Use the model ID enabled in your account.
- Runtime: actual AF3 profile path and the existing library/E034 paths.
- AF3: actual installed Python, runner, model parameters and databases paths. Confirm
  the installation's license prerequisites before setting `license_acknowledged`
  to true. The example paths are placeholders, not an installation procedure.
- Keep `bounded_pair_seeds` false initially. After E037, choose a validated scenario
  and matching chunk settings. Any runtime change requires fresh execution plans.

```bash
read -rsp 'LLM API Key: ' AIDD_LLM_API_KEY
printf '\n'
export AIDD_LLM_API_KEY
python -m json.tool "$AIDD_LLM_PROFILE" >/dev/null
python -m json.tool "$AIDD_RUNTIME_PROFILE" >/dev/null
python -c 'import json,os; from pathlib import Path; from aidd_agent.prediction import load_model_profile; r=json.load(open(os.environ["AIDD_RUNTIME_PROFILE"])); load_model_profile(Path(r["af3_profile"])); print("AF3 profile accepted")'
nvidia-smi
```

Profile acceptance and visible GPUs do not establish that AF3 can execute. The real
prediction below tests the separate AF3 interpreter/environment. Do not put API keys
in JSON, Git, command arguments or returned logs.

## 4. Validate GPT planning, then public protein/structure APIs

First test only live model planning:

```bash
bash scripts/run_e038_prompt.sh --plan-only --prompt 'Find human WEE1 in UniProt, list PDB candidates at resolution at most 3 Angstrom, download PDB 8BJU, and prepare full-length AF3 input with seed 1. Do not run prediction or 3D search.'
```

Inspect the printed `plan.json`: English summary, correct actions and preceding
protein dependencies, no `af3_run` or `search_3d`. A valid JSON response alone is
insufficient if the plan omits requested actions or has the wrong target.

Then run the same request end to end through the data adapters:

```bash
bash scripts/run_e038_prompt.sh --prompt 'Find human WEE1 in UniProt, list PDB candidates at resolution at most 3 Angstrom, download PDB 8BJU, and prepare full-length AF3 input with seed 1. Do not run prediction or 3D search.'
```

Acceptance: report status complete; human WEE1 identity (P30291 / taxon 9606),
FASTA and sequence provenance present; PDB candidate evidence and downloaded mmCIF
present; AF3 input has the verified full-length sequence and seed 1. Prepared input
is not a completed prediction. PDB retrieval is not receptor-quality acceptance.
Each wrapper invocation creates a new plan, including this second call.

## 5. Validate one automatic compute workflow

After steps 1-4 pass, a single real prompt can test preparation, installed AF3 and
both calibrated searches without separately repeating each expensive component:

```bash
bash scripts/run_e038_prompt.sh --allow-compute --prompt 'Retrieve the full-length human WEE1 sequence from UniProt accession P30291, list its PDB candidates at resolution at most 3 Angstrom, download PDB 8BJU, prepare AF3 input with seed 1, run the installed AF3 prediction, then run both calibrated WEE1 queries wee1_qt9 and wee1_824 against the existing full library with Gaussian refinement and E031 annotations only. Use the existing calibrated queries independently of the AF3 prediction; do not use the predicted structure to build a new query.'
```

Acceptance requires ALL requested actions in the plan and ALL steps complete:

- Protein and structure evidence is real and matches the requested identity.
- `af3_run` returns `prediction_completed`, with model, confidence and a hashed
  `prediction-manifest.json`. Input preparation alone does not count.
- `search_3d` covers both calibrated queries (one `wee1_both` action or two individual
  query actions), returns `search_completed`, and child equivalence checks pass.
- Report `e031_changes_ranking` is false; no silent biological-quality claim.

This tests orchestration across supported adapters. The AF3 result does not become
an automatic new search query. New targets, arbitrary receptor-to-query conversion,
multi-protein/NA complexes and unrestricted task execution are not supported by
this prompt interface. Full-length AF3 can require substantial GPU time/memory;
there is no measured end-to-end time estimate for your installation yet.

## 6. Verify completed-stage reuse

Keep the printed context (`db`, `user_id`, `project_id`) and `plan` path from step 5.
Replace the placeholders below with those exact values:

```bash
python -m aidd_agent.prompt_workflow run \
  --db /ACTUAL/registry/aidd.sqlite3 --user USR-ACTUAL --project PRJ-ACTUAL \
  --plan /ACTUAL/Project/runs/PROMPT-ACTUAL/plan.json \
  --runtime "$AIDD_RUNTIME_PROFILE" --allow-compute
```

Acceptance: complete again, original valid artifacts/receipts retained without a
second AF3/search computation. Use this same command to retry an interrupted plan.
Do not invoke the wrapper to resume: it creates a new plan. Do not edit sealed plan
JSON or change code/runtime mid-run; fix configuration before creating a fresh plan.

## 7. Verify boundaries without expensive computation

```bash
# Expected blocked at search_3d, exit code 2, without launching search.
bash scripts/run_e038_prompt.sh --prompt 'Run the calibrated WEE1 QT9 query against the existing full library. Keep E031 annotation-only.'

# Expected clarification about species; no guessed sequence or AF3 execution.
bash scripts/run_e038_prompt.sh --plan-only --prompt 'Prepare AF3 input for WEE1. Ask me for the organism before retrieving a sequence.'

# Expected unsupported-query clarification, never a substituted WEE1 search.
bash scripts/run_e038_prompt.sh --plan-only --prompt 'Run a calibrated 3D search for an EGFR query, not a WEE1 query.'
```

Inspect semantic intent as well as schema validity: a model may produce structurally
valid but inappropriate plans. Treat unexpected behavior as a failed acceptance case.

## Evidence to return

Return the E037 root `report.md`/`report.json`, E038 smoke report, and the live task
`plan.json` plus `execution/report.json`. Include the failed step log if needed,
the AF3 prediction manifest and underlying search report for successful compute.
Record `git rev-parse HEAD`, the chosen model ID, installed AF3 version and hardware.
Exclude secrets and sensitive prompt content. For speed issues also return the
E037 first-chunk profiles. E035 discordant pose exports still need human review.

Completion establishes the supported WEE1 workflow's engineering integration on
this workstation. It does not establish activity enrichment, docking/pose quality,
universal target support or a production latency SLA.
