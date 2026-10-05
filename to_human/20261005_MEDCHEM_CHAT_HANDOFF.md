# MEDCHEM Agent Chat checkpoint - 2026-10-05

Visual follow-up: see 20261005_MEDCHEM_VISUAL_REFRESH.md for the newer two-phase
masthead, cyclic-peptide D and hidden results panel. Its cumulative E099 package
supersedes the E098 package for new installations; the workflow instructions below
remain applicable.

Repository: D:/agent/projects/aidd_structure_guided_chat.
Branch: feature/structure-guided-chat. Git delivery is now user-authorized; the
workstation pull and restart remain user-operated. Continue the existing software
and Chat storage. See the Git delivery section below.

## Actual workstation receipt supplied by the user

The final report is:

```text
/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/plants-parallel20/report.json
kind: block_plants_run
status: complete
jobs: 8445
attempted_jobs: 8445
failed_jobs: 0
scored: 821619
first_job_gate: passed
```

These user-pasted aggregate fields match the expected 273,873 conformers times
three receptors (1RV1_B, 7NA2_A, 5J7F_A). They supersede the October 2 running
checkpoint. Remote files, hashes and poses have not been independently read here.
The earlier pasted `block_plants_prepare / prepared_not_docked` report was preparation
evidence, not the completion receipt. Its ligand_chemistry_reviewed flag was false;
successful docking does not resolve that chemistry review.

Do not resample, repartition, redock, edit frozen signatures or overwrite the outputs.
The expected three partition populations remain E094=384, E095=1188, E096=1189.
N-E rescoring, independent recall, reference cross-docking and affinity validation
remain separate pending work. No rejection of production blocks is authorized by
these exploratory docking scores.

## Implemented in this update

- MEDCHEM Agent title, brand and favicon using D:/agent/med.png unchanged.
- Teal masthead with three locally rendered explicit cyclic-peptide structures:
  pentapeptide, hexapeptide and N-methyl cyclic peptide. They are illustrations,
  not library hits; exact SMILES are in web/peptides.json. Regenerate with
  scripts/build_chat_molecules.py and RDKit; no runtime RDKit dependency for the UI.
- Hidden Workflow coverage opens on hover, keyboard focus or click, closes on
  Escape/outside click, and includes completed CLI docking adoption.
- Prompt example buttons fill the composer without automatically sending.
- Existing Chat sessions, model selector, task queue and PyMOL controls remain.
- New session-owned block_results tool supports attach/list/analyze from ordinary
  prompts, plus an explicit report-attachment form and deterministic slash commands.
- Attachment records the exact source report digest and counters. It does not claim
  output verification or launch an engine. Reports must be inside the active Project;
  the model can attach only a path explicitly present in the current user message.
- Explicit analysis creates a sealed local coordinator plan and runs in the existing
  background queue. It checks the report, output seals, preparation/sample bindings,
  score count and CLI output lock, then calls the existing block statistics adapter.
  The frozen block_plants.py implementation is unchanged.
- Repeated attach/analyze reuses identifiers. Changed inputs and cross-session IDs
  fail. Partial panels remain partial. Existing Cancel/Resume handles analysis jobs;
  it never resumes the external docking engine.
- Analysis reports include absolute artifact paths and ligand chemistry review
  status. Existing task_report/read_artifact tools can inspect the small block summary.
  Large per-pose CSVs remain local artifacts; no N-E engine was added.

## Continue on the workstation

Synchronize the supplied update into the same Group-Agent checkout after stopping
the Chat server and checking that it has no active scientific task. Do not stop or
restart independent completed docking jobs. Do not use a new storage directory.
The update package installer checks old/new file hashes, refuses unrelated local
edits and backs up replaced source files. It does not modify data or start services.

After transferring MEDCHEM_E098_UPDATE.zip to the workstation, extract it into a
separate update directory and run the installer there:

```bash
python install_e098.py --repo /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
python install_e098.py --repo /mnt/medchem_taltio/wrk/yu_agent/Group-Agent --apply
```

The first command only checks. A local-edit/version mismatch stops before replacing
any files. Keep the reported source backup until workstation acceptance is complete.

Restart the existing Chat with its usual environment and credentials:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
bash scripts/run_chat_agent.sh \
  --storage-root /mnt/local/hand/yuzhang/aidd/e097-chat-workspace \
  --runtime /mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/e097-config/runtime.json \
  --llm-config-dir /mnt/local/hand/yuzhang/aidd/config/llm \
  --port 8766 --allow-compute
```

Use the printed access link and existing conversation. `--allow-compute` retains the
previous workflow capabilities; analysis of attached scores itself does not require
PLANTS, SPORES, a license check, a model API call or a new full-library computation.
Old code-bound unfinished workflow plans must not be blindly resumed after updating.

First prompt (copy the complete path):

```text
Attach this completed PLANTS docking report to this conversation. Show its status,
job count, failed jobs and scored conformer/receptor pairs. Do not rerun docking:
/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/plants-parallel20/report.json
```

Then:

```text
Analyze the attached completed docking result. Compare E094, E095 and E096 separately
for each receptor, including score distributions, missing scores and population
weighting. Preserve CID and pose provenance. Do not resample or redock.
```

After the queued task completes:

```text
Read the completed block analysis and summarize the differences between the three
partition schemes for each receptor. Explain what the evidence supports and what
remains unvalidated. Show the summary CSV and pose handoff paths.
```

Provider-independent commands are also available:

```text
/dock_attach /absolute/path/to/plants-parallel20/report.json
/dock_results
/dock_analyze ATTACHMENT_ID
/status TASK_ID
/results TASK_ID
```

The attachment ID is returned by the first command, not the original preparation
task ID. Use actual IDs, never copy placeholders literally. No final block-analysis
results from the real 821,619-row panel have been received yet.

## Validation and next work

Protocol: experiments/E098-medchem-chat.md. Independent synthetic fixtures check
adoption, duplicate dispatch, Project/session ownership, input mutation, active-run
locking, natural-language tool routing and the real background Python executor.
Browser evidence: data/e098-ui/report.json and desktop.png, coverage.png, mobile.png.
Desktop 1440x1100 and mobile 390x844 checked; actual browser rendering inspected.
Live LLM routing and actual 821,619-row analysis acceptance remain workstation checks.

Full local regression: 743 passed, 6 skipped (optional dependencies), 75.27 seconds.
The minimal virtual environment lacked requests; the final run reused the existing
local requests/urllib3/idna/certifi/charset_normalizer packages in an isolated ignored
test dependency directory. It did not change the workstation environment or scientific
implementations. An earlier broad conda fallback exposed an incompatible optional
Numba; the final isolated dependency setup avoids mixing that backend. English-content
and JavaScript syntax checks passed. No live provider or external engine ran in QA.

Next: synchronize/restart Chat, attach the final report, run and inspect real block
analysis, review chemistry and nearby removed receptor components, then define the
next evaluation/N-E stage from the actual outputs. Do not infer improved retrieval
recall or affinity from successful PLANTS execution.


## Git delivery after UI approval

The user approved the MEDCHEM design and Helvetica typography and requested remote
synchronization. Pull the existing feature branch after stopping the Chat server:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
```

Use the restart command above, retain the e097-chat-workspace storage and open the
printed access link. Do not install the ZIP package on top of this Git update.
No scientific engine rerun is part of this update. First attach plants-parallel20's
completed report and inspect 8445 jobs / 0 failures / 821619 scores. Then request
per-receptor E094/E095/E096 analysis. Keep the analysis task ID and inspect its output;
there may be substantial I/O while checking all original file hashes.
