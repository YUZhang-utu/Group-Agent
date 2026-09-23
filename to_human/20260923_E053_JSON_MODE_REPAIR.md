# E053 recommendation JSON-mode repair - 2026-09-23

The user reported `consensus_recommend` failing with HTTP 400,
`invalid_request_error`, 48031 prompt characters and a 120-second timeout.
The report's `local-coordinator` model identifies the local workflow coordinator,
not the remote provider/model. The error alone does not establish context overflow.

The consensus-specific system prompt did not explicitly request JSON, while the
shared client enabled `response_format: {"type": "json_object"}` by default.
[DeepSeek's JSON-mode contract](https://api-docs.deepseek.com/guides/json_mode/)
requires the word `json` in a system or user message. The shared client now adds
an explicit single-JSON-object instruction for every caller, including custom
consensus prompts and profiles with provider JSON mode disabled. Local strict
JSON parsing and scientific design validation remain required.

Regression protocol: use a mock provider that rejects JSON-mode requests without
the required instruction, then exercise the actual recommendation request with
126 anchors and 42 templates. Before the fix, it reproduced HTTP 400 with
`invalid_request_error`; after the fix, the recommendation succeeds with complete
evidence. Two additional cases cover custom prompts with JSON mode on/off and
rejection of invalid output by the local validator. Focused validation:

```powershell
$env:PYTHONPATH='src'
python -m pytest tests/test_prompt_workflow.py tests/test_consensus.py tests/test_container_providers.py -q
python scripts/check_english.py
```

Result: 41 passed, one skipped; English-content guard passed. This confirms the
request compatibility repair locally, not the cause or resolution of the actual
workstation response. No live provider call or full-library search was performed.
Evidence is not truncated and no provider/model fallback is introduced.

## Workstation continuation after receiving the updated branch

Stop the existing Chat process with Ctrl+C in its server terminal. Preserve the
existing runtime/LLM profile environment variables and Chat storage, then run:

```bash
conda activate aidd-workstation
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
git log -3 --oneline
export AIDD_PY="$(command -v python)"
bash scripts/run_chat_agent.sh \
  --storage-root /mnt/local/hand/yuzhang/aidd/e053-chat-workspace \
  --runtime "${AIDD_RUNTIME_PROFILE:?AIDD_RUNTIME_PROFILE is not set}" \
  --llm-config-dir "${AIDD_LLM_CONFIG_DIR:?AIDD_LLM_CONFIG_DIR is not set}" \
  --port 8766 --allow-compute
```

Stop dependent steps if an update/environment command fails; preserve local edits
and untracked files. In the same Chat conversation, send a new `/recommend`.
Reuse completed discovery/consensus; do not resume the failed code-hashed task.
If it succeeds, review the design before `/adopt`, then `/guided` for the intended
consensus funnel. The previously stopped E050 job remains stopped.

If HTTP 400 persists, capture the new error, task ID and checkout commit, plus
the actual selected provider/model and JSON-mode setting (never the API key).
The generic allowlisted error code cannot distinguish all provider rejection
causes. Do not shorten scientific evidence or change models on that basis alone.

The prior workstation paths and scientific checkpoint remain documented in
[the overnight handoff](20260922_E053_OVERNIGHT_HANDOFF.md).
