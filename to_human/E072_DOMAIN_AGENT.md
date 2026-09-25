# Domain agent: update and workstation acceptance

Ordinary Chat messages now enter a bounded tool-feedback loop. The model can inspect
existing tasks, read report fields, query the configured molecule registry, retrieve
Europe PMC abstracts, invoke existing workflow adapters and use PyMOL. Each result
returns to the model before it decides the next action. Existing slash commands
remain available. This is an integration release, not evidence of superiority to
Codex or guaranteed understanding of every prompt.

## Update

In the workstation repository, check `git status`, preserve any local changes, then:

```bash
git pull --ff-only origin feature/structure-guided-chat
bash scripts/run_chat_agent.sh \
  --storage-root "$AIDD_NEW_WORKSPACE" \
  --runtime "$AIDD_RUNTIME_PROFILE" \
  --llm-config-dir "$AIDD_LLM_CONFIG_DIR" \
  --port 8766 --allow-compute
```

Stop the previous Chat server before launching its replacement. Retain the same
workspace, runtime and provider configuration to retain existing sessions/tasks.
Reload the browser to load the updated JavaScript. Do not restart a completed
full-library search for these checks. Select the session containing the task.

The library connection reuses `search.batch` inside the existing runtime JSON:

```json
{"search": {"batch": "/mnt/local/hand/yuzhang/aidd/library-precompute-20260911"}}
```

Merge this field into the existing configuration; do not replace other settings.
The directory must contain `artifacts/catalog.json` and `registry.sqlite3`.
Scientific search adapters also use the existing chemical companion artifacts.
Provider credentials remain in the existing provider configuration/environment.
Literature retrieval requires outbound HTTPS to the Europe PMC service.

## Natural-language acceptance prompts

Replace TASK_ID and MOL_ID with real identifiers. These are test requests, not
claims that the live provider has already completed them.

1. `Check which molecule library is connected. Report its catalog size and whether the molecule registry and chemical companion catalog exist. Do not start a search.`
2. `Find molecule MOL_ID in my connected library. List its original name, conformer IDs, source MOL2 files and record indices. Keep molecule and conformer identities separate.`
3. `Inspect completed task TASK_ID. Explain the actual screening funnel, retained molecule count and validation limitations. Do not rerun it.`
4. `Read task TASK_ID and search published MDM2 macrocycle research. Compare the reported observations with the retrieved evidence, cite the sources, and distinguish hypotheses from established conclusions.`
5. `Open the aligned reference structures from task TASK_ID in PyMOL. Report whether the viewer operation completed or is still pending.`
6. `In the current PyMOL viewer, inspect the structure IDs and show whole protein residues within 5 angstrom of the ligand in v001. Keep the ligand visible as sticks, color by element, and label the nearby residues.`
7. `Inspect the supported interaction types for v001. Show only geometrically detected interactions, use different colors by type, and explain missing chemistry metadata and unsupported types.`
8. `Explain the AF3 confidence results from task TASK_ID, including ipTM and any available interface metrics. Explain what they cannot tell us about binding affinity.`
9. `List the available workflows and their implementation or validation limitations. Tell me which ones require installed external tools or a connected desktop viewer.`

Existing PDB acquisition, shared-pocket cocrystal collection, alignment, ligand
diversity, consensus design, library screening, budget pages/MOL2 export, AF3
CCD/SMILES preparation and prediction, interaction analysis, and configured docking
remain accessible through the existing workflow decision interface. They retain
their original prerequisites, validation status and compute policies. The new loop
does not turn unavailable adapters into implemented methods.

## Inspect results

The task panel shows agent runs and individual tool steps, purposes, status and
errors. Download the trace from its card. Traces are also stored under
`STORAGE_ROOT/chat/agents/SESSION_ID/RUN_ID.json`. Scientific tasks retain their
existing task IDs, reports, logs and output paths. Queued computation is asynchronous;
ask for its status later. Viewer state and receipts remain in the existing viewer
directory. Read-only analysis must not create a new scientific task.

## Limits and verification

- At most ten planning steps per request; a time check follows tool execution.
  This is not a strict wall-clock timeout on an in-flight provider request.
- Pending dispatch blocks further mutation in that request. Resume returns the
  existing task ID; duplicate operations are blocked within the request.
- Reports are restricted to the session's owned project; artifact reads use
  discovered IDs rather than arbitrary model-selected filesystem paths.
- Registry lookups currently use exact molecule ID or exact original source name.
  They do not provide arbitrary SQL, fuzzy search or automatic MOL2 extraction.
  Use the existing export/budget-page workflow for candidate file export.
- Literature support is title/abstract search, not unrestricted browsing or
  full-text systematic review. Retrieved evidence is not instruction text.
- PyMOL uses the existing bounded API interpreter and scene metadata; the model
  does not inspect screenshot pixels. Existing chemistry limitations still apply.
- Saved agent traces are auditable but the agent loop itself does not automatically
  resume after interruption. Scientific job resume is a separate existing feature.
- Local tests use controlled planners and tool fixtures. Real-provider understanding,
  workstation library availability, literature networking and desktop rendering
  must be checked with the prompts above before a demonstration.

For a quality benchmark, use ten paraphrases of each key request and compare task
completion, correct selections, cited evidence, unexpected jobs, time and cost
against the previous router with identical inputs. Do not infer quality from one
successful prompt or from fixture tests alone.
