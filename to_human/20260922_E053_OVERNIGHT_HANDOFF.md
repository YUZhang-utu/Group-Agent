# E053 overnight handoff - 2026-09-22

## Current checkpoint

The user will verify tomorrow. The workstation reached `consensus_recommend`,
but no successful live recommendation has been reported. The latest reported
failure is only `RuntimeError`; its underlying cause remains unknown. Do not
assume a timeout, invalid credential, quota problem or context-size rejection.

Latest executable fix pushed: **f64bbb8**, branch **feature/structure-guided-chat**.
It preserves safe HTTP/timeout/connection diagnostics in the workflow report.
The user has not yet confirmed pulling or testing this commit. This document is
a later documentation checkpoint, so the branch HEAD may be newer than f64bbb8.
No new E053 full-library scan is confirmed running; there is no overnight search
result to expect yet. The earlier E050 run was explicitly stopped by the user.

## Workstation location and environment

- Main working checkout: `/mnt/medchem_taltio/wrk/yu_agent/Group-Agent`.
- Branch: `feature/structure-guided-chat`.
- Old `/mnt/medchem_taltio/wrk/yu_agent/Group-Agent-e051` worktree was detached
  successfully to release the branch. Keep it; there is no need to delete it.
- The main checkout had an untracked `k`; its contents/type are unknown and it
  was preserved. Do not run `git clean` or `reset --hard` to remove it.
- Scientific environment: `aidd-workstation`. A later shell showed `base`;
  this alone does not prove which environment the running Chat used.
- Missing dependency `requests` was reported; installation in the scientific
  environment was instructed. Later progress reached recommendation, but there
  is no separate pasted installation receipt.
- Chat storage: `/mnt/local/hand/yuzhang/aidd/e053-chat-workspace`.
- Project observed in a failed diversity task:
  `prj-2a23d9a15368-prompt-aidd` (user `workstation`).
- Port: `8766`. Runtime and LLM profiles are supplied by the existing
  `AIDD_RUNTIME_PROFILE` and `AIDD_LLM_CONFIG_DIR`; their concrete workstation
  paths have not been supplied here. Do not invent replacement paths or keys.
- Chat workers inherit the Python interpreter that started the Chat service.
  Activating another environment in another terminal does not change that service.

## Target and intended workflow

This is a general-target implementation tested with human MDM2, Q00987 / 9606.
Reference: `5C5A:NUT:A:201`, target protein author chain `A`.

The user wants all scientific interaction through Chat, with human review of
outputs. The current implementation uses separate tasks, not a verified automatic
overnight chain from one long prompt:

`structure_diversity` -> `structure_consensus` -> `/recommend` -> `/adopt` -> `/guided`

Wait for each task to complete. The final `/guided` launches the full configured
catalog, not a library pilot. Human anchor-combination selection follows search.
Do not accidentally route this to WEE1, the stopped E050 run, or benchmark_funnel.
The fact that the user reached consensus_recommend indicates a completed source
consensus was accepted by the coordinator; its successful workstation run ID and
full report have not been pasted. Reuse that source in the same conversation.

## Failure history and fixes

| Observation | Diagnosis / status | Commit |
|---|---|---|
| Branch already used by Group-Agent-e051 | Detached old worktree; switched and updated main Group-Agent successfully | ee58141 confirmed locally by user |
| benchmark_funnel: coarse audit needs explicit constraints | Separate wrong-path task; not evidence of an MDM2 diversity failure | No need to rerun this audit |
| protein_fetch complete; structure_diversity ModuleNotFoundError | Missing module was hidden; user subsequently identified requests | d9f7ba5 reports missing module and worker interpreter |
| consensus_recommend: prompt must contain 1-20000 characters | Consensus allowed 100000 but common request helper allowed only 20000 | 3445faf aligns the trusted evidence budget without removing evidence |
| consensus_recommend: RuntimeError | Workflow hid safe transport details; exact provider fault unknown | f64bbb8 preserves HTTP status, selected provider codes, timeout and connection diagnostics |

Known failed diversity report (historical, not the current recommendation):

`/mnt/local/hand/yuzhang/aidd/e053-chat-workspace/users/workstation/projects/prj-2a23d9a15368-prompt-aidd/runs/PROMPT-d5a3e840acfd4196/execution/report.json`

The latest recommendation task ID/report path has not been provided. Obtain it
from the same Chat task list rather than guessing or reusing the historical path.

## Tomorrow: exact continuation

1. Stop the old Chat server with Ctrl+C in its terminal, if it is still running.
2. In that terminal, update the main checkout and confirm the interpreter:

```bash
conda activate aidd-workstation
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
git log -3 --oneline
export AIDD_PY="$(command -v python)"
"$AIDD_PY" -c 'import sys; print(sys.executable); import rdkit, gemmi, Bio, numpy, scipy, requests; from rdkit.Chem import rdFingerprintGenerator; print("Dependencies OK")'
```

If a command fails, retain the output and stop the dependent steps. Do not erase
local edits, outputs or profiles. Do not reinstall the whole environment.

3. Restart using the same storage and existing profiles:

```bash
bash scripts/run_chat_agent.sh \
  --storage-root /mnt/local/hand/yuzhang/aidd/e053-chat-workspace \
  --runtime "${AIDD_RUNTIME_PROFILE:?AIDD_RUNTIME_PROFILE is not set}" \
  --llm-config-dir "${AIDD_LLM_CONFIG_DIR:?AIDD_LLM_CONFIG_DIR is not set}" \
  --port 8766 --allow-compute
```

4. Open the **same conversation**, and send `/recommend` to create a new task.
   Do not resume a failed code-hashed task after updating code. Do not redo the
   completed structural discovery or consensus.
5. If it fails, copy the new task's complete `error`, task ID, report path and
   the checkout commit. No API key or credential-bearing profile contents are
   needed. If the error is still only RuntimeError, check the actual server
   checkout/interpreter and investigate other runtime sources rather than
   assuming the provider failed.
6. If recommendation succeeds, review it, send `/adopt`, wait for completion,
   then `/guided`. Confirm the action is `consensus_funnel` and actual catalog
   progress is advancing. Starting Chat or completing recommendation alone does
   not launch the full-library search. Keep the service running during execution.

## Evidence already available, and its limits

- Implementation baseline: 372 tests passed, two skipped.
- Latest diagnostic/context-focused suite: 34 tests passed, plus English-content
  guard. These are local tests, not a successful live-provider recommendation.
- Local real MDM2 preparation: 48 admitted instances from 45 PDBs; 42 distinct
  prepared ligand templates; 126 spatial/directional interaction modes; one met
  the automatic mandatory-support floor; eight templates in the test design.
- Deterministic injected recommendation control: 43/48 crystal identity poses
  passed. This is not a live LLM recommendation or independent activity validation.
- Independent binding panel: 24 ChEMBL molecules, one reference-graph overlap
  excluded, 23/23 remaining controls retained in the exploratory funnel.
- Detailed local evidence: `data/e053-mdm2-implementation.json` and
  `D:/agent/MDM2/analysis/e053` on the development machine, not the Linux workstation.
- Full-library coverage, rejection fractions, latency, matched-negative enrichment
  and general-target recall remain unverified on the workstation.
- Polymer chemistry, transformed assemblies, other receptor states, automatic
  throughput optimization, and general-target docking/affinity adapters are not
  completed capabilities. The current supported path is nonpolymer references
  in a conservatively admitted receptor-state cohort.

## Continuity instruction

Read this handoff, research-state.yaml and the latest research_log.md entries.
Continue from the failed live consensus recommendation, not from PDB discovery.
First use the new diagnostic to identify the actual provider/runtime failure.
Preserve the user's main checkout, existing Chat conversation, completed artifacts
and untracked k. Do not launch duplicate scans or claim the full-library task ran.
