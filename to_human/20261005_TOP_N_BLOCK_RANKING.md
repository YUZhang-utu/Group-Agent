# Top-N block ranking and candidate inspection

User-confirmed default: ten DIFFERENT molecules per block, using each molecule's
best in-block conformer for a given receptor, then averaging those ten scores.
Lower ChemPLP is better. Whole-block means remain diagnostics, not ranking keys.
Top-5 means are also available. Custom N from 1 to 100 is supported.

Local validation: 768 tests passed, 6 skipped. Pull the existing
`feature/structure-guided-chat` branch and restart Chat in the same workspace.
This delivery also includes the E101 reserved rescoring contract, not an installed
N-E, GNINA or affinity execution adapter.

## In the existing Chat Project

Attach the completed PLANTS report once, as described in the MEDCHEM handoff.
Then request:

> Analyze the attached completed docking result. Rank blocks by the mean of their
> best ten distinct molecules, using each molecule's best in-block conformer.
> Show the five leading blocks separately for each partition scheme and receptor.

After the queued task finishes:

> Show the Top-5 molecules of the leading E096 block for receptor 5J7F_A, including
> molecule IDs, conformer IDs, scores and saved pose references.

Use actual group/block identifiers returned by the tool; never invent IDs.
Existing pre-update analysis tasks remain historical. Request a new analysis to
obtain Top-N outputs; the docking engine is not invoked again.

Provider-independent commands:

```text
/dock_results
/dock_analyze ATTACHMENT_ID 10 molecule
/block_ranks TASK_ID
/block_ranks TASK_ID E096 5J7F_A 10
/block_ranks TASK_ID E096 5J7F_A 5
/block_top TASK_ID E096 5J7F_A BLOCK_ID 5
```

The first ranks query discovers the actual scheme/receptor groups. A changed N
or ranking unit creates a separate analysis task. Repeated identical requests
reuse a queued/running/completed task. Failed analyses may be queued afresh.

## Outputs

- `block_rankings.csv`: Top-5, Top-10 and requested-N means, within-group ranks,
  molecule/conformer counts, missingness and eligibility.
- `block_top_candidates.csv`: retained leading distinct molecules, selected CID,
  score, pose name, file and original pose record index.
- `block_ranking.sqlite`: indexed bounded Chat lookup; verified against its seal.
- `ranking_review.md`: primary scoring rule and leading blocks.
- `block_summary.csv`, `block_scores.csv`, `review.md`: previous diagnostics.

Candidates are separate within each block/receptor. The same molecule may occur
in several blocks; no independence or additive hit-count claim follows. All
scores currently come from the existing ChemPLP run. N-E code is pending and no
GNINA or affinity inference is launched by these commands.

Blocks with fewer than N successful units have a null Top-N mean and no primary
rank; they are not padded or silently averaged over fewer units. Incomplete
panels are excluded from the primary ranking. Their saved candidates can still
be inspected with explicit status. Sampling budgets differ between schemes;
ranking useful blocks does not itself prove one partition scheme superior.

## CLI alternative

```bash
python -m aidd_agent.block_ranking \
  --source /absolute/path/plants-parallel20/report.json \
  --output /absolute/path/new-top10-analysis \
  --top-n 10 --ranking-unit molecule
```

Use a fresh output directory and the existing project Python environment. The
CLI verifies saved docking output seals and does not execute a model. Large
pose collections may take time to hash. Remote result values have not been
inspected locally; the workstation must execute this analysis or the Chat task
before any actual winning block IDs can be reported.

## Recovery from a planner argument error

Before restarting, use `/status` to check the current task; let queued/running
work finish rather than interrupting it. After updating, restart the Chat server
process, not only the browser. Run
`/dock_results` first. The response now includes `analysis_tasks` with IDs/status
and `argument_contract_version: top-n-arguments-v1`. An old response without this
version indicates the server has not loaded this recovery revision.

- If an analysis task is queued/running, inspect `/status TASK_ID`; do not submit
  another analysis. Completed tasks can be inspected with `/block_ranks TASK_ID`.
- If only an attachment exists, use `/dock_analyze ATTACHMENT_ID 10 molecule`.
- If neither exists, attach the exact original final PLANTS report with
  `/dock_attach /absolute/path/to/report.json`, then analyze its returned ID.

The planner receives machine-readable operation-specific field contracts.
Unknown fields are rejected with their exact names and allowed alternatives;
they are not silently ignored or translated. Successfully queued block analysis
returns its task receipt without another planner turn. Identical repeated tool
failures stop after two occurrences while retaining the trace and prior receipts.
The actual workstation error arguments require its Agent trace; local tests
demonstrate recovery behavior but do not establish which field failed remotely.
