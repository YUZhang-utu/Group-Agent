# Saved scores to block-restricted search in MEDCHEM Agent

Implementation: E104. Local fixture validation; real workstation search acceptance
is pending. Existing full EquiScore/PLANTS outputs remain unchanged.

## Update and open

In the workstation repository:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
```

Restart the existing Chat server with its original workspace, user, Project,
runtime and Python environment. Keep `--allow-compute` for search/docking. Do not
start a new workspace or reinstall the working EquiScore environment. Refresh
the browser after restart. Use the original MDM2 conversation so its query task
remains available.

## Import the existing result (no repeated inference)

Open **Tasks & results → Attach an existing result**. Choose **Completed EquiScore
analysis**. Supply the final analysis report:

```text
/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/equiscore-v3/equiscore-analysis/report.json
```

Use the actual path if the analysis was saved elsewhere. The importer requires
the analysis, full EquiScore run and original PLANTS report inside this Project.
It verifies output hashes, full-run success, and the matching ChemPLP baseline.
Import is a queued task because verifying large saved tables can take time.
It copies compact inspection artifacts and keeps links/hashes to the original
panel. It does not run EquiScore, PLANTS or a new block analysis.

Equivalent chat prompt:

> Import my completed full EquiScore analysis from [paste the exact report path].
> Record the existing rankings and ChemPLP comparison without repeating inference.

## Select blocks and launch

In **Tasks & results → Leading blocks to 3D search**:

1. **Load saved analyses and queries** after import completes.
2. Choose the imported analysis and the existing confirmed MDM2 query task.
   If no query appears, use its original conversation or attach its existing
   sealed Chat plan. A receptor-assessment or sampling task alone is not a 3D query.
3. Choose one partition and receptor. Select EquiScore, ChemPLP, their union or
   intersection, with Top 5 or Top 10 blocks per score.
4. **Preview leading blocks** to see both rankings, selected IDs and total block
   population. Scores are never pooled across receptors or partition schemes.
5. Set **100000 conformers**, then **Search and dock conformers**.

The budget covers the selected block union, not each block. Multiple conformers
of a molecule are permitted. A source conformer shared by scoring methods is
searched/docked once. The block ranking statistic still uses the mean of the
best 10 distinct molecules; the search-delivery counting unit is different.

Equivalent prompt (replace task IDs with actual returned IDs):

> Using imported analysis task ANALYSIS_ID and confirmed MDM2 query task QUERY_ID,
> select the union of the Top-5 EquiScore and Top-5 ChemPLP blocks for E095 / 7NA2_A.
> Search only those blocks and retain up to 100000 conformers, allowing multiple
> conformers per molecule. Then dock them against 7NA2_A and record the results.

Deterministic fallback, without a model call:

```text
/block_campaign {"operation":"list"}
/block_campaign {"operation":"preview","task_id":"ANALYSIS_ID","scheme":"E095","receptor":"7NA2_A","method":"union","blocks":5}
/block_campaign {"operation":"start","task_id":"ANALYSIS_ID","query_task_id":"QUERY_ID","scheme":"E095","receptor":"7NA2_A","method":"union","blocks":5,"conformers":100000,"dock":true}
```

The IDs must be actual 16-character task IDs, not run IDs or the literal examples.
Use `/status TASK_ID` and `/results TASK_ID` for progress and saved artifacts.
Repeated identical requests return the existing task; failed tasks are not
silently resubmitted. Inspect the log before explicitly resuming.

## What the continuation actually does

- Binds selected memberships to the original sample's frozen database receipts.
- Maps source paths/record indices and content hashes to search artifact IDs.
  Missing or ambiguous mapping stops the run; no outside-block fallback.
- Uses the existing contact-first 3D engine and selected MDM2 query. It refines
  the selected block population, which can exceed the 100000 delivery budget.
  Per-template contact ranks (Gaussian tie-break) are fused by conformer RRF.
- Retains at most the requested number of qualified source conformers. It reports
  shortfall explicitly rather than padding with outside-block entries.
- Exports original MOL2 atom/bond chemistry and records search transforms.
  PLANTS performs its own pose search; input coordinates are the original source
  coordinates, not a claim that the search transform is a docked pose.
- Inherits the original reviewed receptor, binding site, ligand mode and PLANTS
  settings. The trusted configured executable must match the original engine
  hash; only scheduling/timeouts come from the current runtime profile.
- Saves `docking_candidates.csv` with source conformer/molecule/block IDs,
  ChemPLP scores and saved pose references, and Top-10 conformers in the report.

One selected receptor and 100000 delivered conformers produce 100000 planned
conformer/receptor pairs. Distinct-molecule count is reported separately.
Other scheme/receptor groups can be submitted as separate campaigns.

Required existing runtime entries: `search.batch`, `search.workers`,
`search.refine_chunk`, `block_evaluation.sampling_profile` and `plants_profile`.
No EquiScore environment change is needed. New outputs use a fresh sealed Chat
run, so successful old code-hashed calculations remain intact.

These are exploratory priorities. Follow-up panels are search-selected, not
uniform block samples; no population estimate, measured affinity or experimental
enrichment is implied. This continuation does not automatically run a second
EquiScore inference on the newly docked candidates.
