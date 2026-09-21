# E041: conversational AIDD workbench

The browser is the entry point; existing scientific scripts run in the background.
It supports persistent conversations, GPT/DeepSeek selection, background tasks,
actual report paths, cancellation and explicit resume. Local UI and maintained
responses are English; user requests can be in any language the selected model
understands. Provider planning remains subject to semantic validation.

## Start once on the workstation

Activate the existing AIDD environment and keep the configured API keys in your
shell environment. Use the same local provider/runtime profiles as E040:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
export AIDD_LLM_CONFIG_DIR=/mnt/local/hand/yuzhang/aidd/config/llm
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime-apptainer.local.json
bash scripts/run_chat_agent.sh --allow-compute
```

Open the printed `http://127.0.0.1:8765/#token=...` link in the workstation browser.
The local access token is not an LLM API key. It is required for the local API,
removed from the visible URL after loading and kept in browser session storage.
Keep the server running (for example inside your existing tmux session). Closing
the browser does not cancel tasks. Closing the server interrupts active execution.
Only one server may use this storage directory at a time.

For access from another computer, use your normal SSH connection with
`-L 8765:127.0.0.1:8765`, then open the same local link. The server intentionally
binds to loopback; it is a single-operator workstation tool, not a public multi-user
service. Do not expose it directly to the internet.

Without `--allow-compute`, retrieval/preparation remain available, but scientific
compute reports blocked. There is no new per-task approval dialog after starting
with compute enabled. Select DeepSeek or GPT in the page; only that provider's key
is required in the server environment.

## Work in one conversation

1. Create a conversation and ask: `Retrieve human WEE1 P30291 and predict its
   full-length structure with AF3, seed 1.`
2. Ask: `What is the status of my task?` while it runs. The task panel also updates.
3. After completion ask: `Where are the model and report files?`
4. Ask: `Run both calibrated WEE1 searches and keep E031 annotation-only.`
5. Ask: `Where is my latest search report?`

Status/results answers are assembled from actual task records/reports; they do not
start another scientific job. The task panel displays the plan, report and log paths
with copy buttons, step records and the execution log. Progress is stage/log based,
not an invented percentage or completion-time estimate. Completed tasks append a
message automatically. Heavy jobs run serially to avoid accidental GPU/CPU overload.

Missing species or unsupported methods should trigger an English clarification.
Your reply remains in the same conversation. Full history is persisted; the model
receives a bounded recent conversation and task summaries, not all historical files.
The model receives more context than the old one-shot planner: recent conversation
text may contain paths and result summaries. Raw libraries, structure files and
full execution logs are not automatically uploaded by the router.

Deterministic commands are available even without a working model connection:

```text
/status
/results
/capabilities
/cancel
/resume
```

They use the latest task, or accept a task ID after the command. Buttons use the
same controls. Cancel stops the active process group and retains artifacts. A
cancel during model planning waits for the current HTTP request to finish but
prevents launching compute. Resume queues the same sealed plan and receipt checks;
changed code/configuration still requires a new task. Resume a clarification-only
plan by answering its question in a new message instead of replaying that plan.
Interrupted tasks are not automatically relaunched after server restart.

## Bring existing results into the conversation

Use `Attach an existing result` in the task panel. Paste the owned run's `plan.json`
path. No scientific computation is launched. Plans must belong to this active
Project and pass the stored seal; completed reports must match their receipt.
This checks report provenance, not independent scientific validation of every file.

For your existing AF3 result:

```text
/mnt/local/hand/yuzhang/aidd/prompt-workspace/users/workstation/projects/prj-252fa94197d6-prompt-aidd/runs/PROMPT-184f88ad866840d7/plan.json
```

For your existing full-library search:

```text
/mnt/local/hand/yuzhang/aidd/prompt-workspace/users/workstation/projects/prj-252fa94197d6-prompt-aidd/runs/PROMPT-4689d49ad1834f89/plan.json
```

Then ask where the results are, or use `/results`. Old runs may not be resumable
after a code update; attachment is read-only and does not require rerunning them.

## Scientific scope and what remains

| Stage | Current implementation | Acceptance still required |
|---|---|---|
| Protein/PDB evidence | Generic verified accessions/species and PDB retrieval | Task-specific identity and receptor choice |
| AF3 | Generic single-protein/CCD input and installed container execution | Structure/region pLDDT/PAE and construct review; confidence is not activity |
| 3D retrieval | Full-library USRCAT search for two calibrated WEE1 queries | New-query calibration and independent held-out panels |
| Gaussian refinement | Rigid shape/feature refinement of candidate conformers | Pose-quality validation; no torsion relaxation claim |
| E031 | Equivalent interaction annotations | Human pose review; remains annotation-only |
| Molecule aggregation / pocket QC | Existing CLI components | Chat adapters with explicit site/query artifact bindings |
| Glide / PLANTS docking | User reports installed tools; no validated executor in chat yet | Executable/license discovery, grid/receptor/ligand prep, redocking and cross-docking |
| Rescoring / diversity / properties | Not a completed general pipeline | Chosen methods, reference datasets and calibrated selection criteria |
| Biological validation | Not established | Held-out active/decoy labels and ultimately experimental evidence |

Do not interpret the chat interface as completing missing scientific stages. It
never turns an arbitrary AF3 protein into a calibrated ligand query. User-reported
current search execution: QT9 32.859 s; 824 64.056 s, fresh computation and exact
Gaussian/E031 equivalence. Startup/index load and final checks are separate; the
E037 controlled speed comparison remains pending.

## Discover installed docking tools without running a job

After your usual Schr?dinger/Maestro module load, run:

```bash
bash scripts/probe_docking_tools.sh
```

If the module command is available to the script, you can instead supply its exact
site module name as the sole argument. This lists the SCHRODINGER root and available
Glide/LigPrep/PrepWizard/PLANTS paths; it does not launch docking or test licenses.
Return the JSON output. Maestro being available does not itself establish Glide
execution or a validated docking protocol. Next integrate trusted engine profiles,
then receptor/grid preparation and a fixed reference-ligand redocking protocol,
then candidate docking with provenance and separate timing/quality gates.

## Workstation chat acceptance

Evidence classification, explicit same-pose selection previews and a separate
export confirmation are now available through the [E042 workflow](E042_SCREENING_TO_DOCKING_HANDOFF.md).
These reuse completed searches and do not submit docking jobs. Compact anchor
evidence summaries are included in model context for selection follow-ups.

- Refresh and reopen: sessions, messages and task paths remain available.
- Ask status/results repeatedly: no duplicate compute job is created.
- Run one existing calibrated search; verify both inner equivalence checks.
- Reuse the completed task with `/resume` before changing code/configuration:
  existing valid stage receipts should avoid repeating expensive computation.
- Test a missing organism and an unsupported new target: clarification, no WEE1 substitution.
- Cancel a disposable run: process stops, files remain, task is cancelled.
- Attach the two earlier completed runs: results can be located without rerunning.
- Compare GPT and DeepSeek on the same request; preserve plan correctness evidence.

Local tests cover persistent state, ownership, real child-process clarification
execution, cancellation, HTTP authentication/origin handling, report-derived answers
and invalid model decisions. They do not replace real model/GPU/browser acceptance.
