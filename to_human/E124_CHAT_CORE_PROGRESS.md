# AIDD core workflow audit and integration

Updated: 2026-10-29. Historical E124 snapshot; the E125 general workflow delivery
includes these changes and supersedes the execution boundaries listed below.
See `E125_GENERAL_PROMPT_WORKFLOWS.md` for the current interface.

## Current direction

Defer block-routing and block-benchmark development. E109/E110 did not establish
a consistent large routing benefit. Preserve code, experiments and saved results;
historical inspection and explicitly requested block workflows remain available.
Chat now states this priority rather than recommending blocks by default.

The active priorities are AF3, target-specific full-library screening and reliable
Chat access to task state and saved evidence.

## What is connected

| Component | Current capability | Evidence and remaining boundary |
| --- | --- | --- |
| Chat task execution | Persistent sessions, typed plans, background queue, status/results, cancellation and guarded resume | Offline regression coverage; actual worker environment and runtime profiles remain host-specific |
| AF3 preparation and execution | Verified protein, requested construct, seeds, CCD and literal SMILES; native/Apptainer runtime | Chat now supports explicitly requested identical protein copies, 1..8; heteromer source selection remains outside this adapter |
| Independent AF3 ligand batches | Standalone batch runner shares protein MSA/templates and runs ligand inference sequentially | Batch execution remains CLI-only; report evidence can now enter Chat |
| Full-library Chat workflows | Existing `full_count`, `funnel`, guided/consensus, selection and export adapters | These have different search scopes. ANN budgets and descriptor refinement budgets must not be called uncapped pose evaluation |
| KRAS production workflow | Completed E123 pool scoring and MOL2 export through standalone scripts | Production orchestration is still CLI-only. Chat can now inspect the final report without rerunning the search |
| External result bridge | KRAS E123 report, AF3 batch report and AF3 prediction manifest import; paginated reading and explicit source freshness check | Report schema/counter validation only; remote MOL2/CIF payloads are not opened or hash-verified by import |

## Verified KRAS progress

The returned final archive and individual receipts were audited separately. The
new Chat bridge was also exercised against the actual returned final report.

| Target | Existing pool IDs processed | Balanced pool | Selected IDs | Exported MOL2 records | Quarantined IDs | MOL2 files |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KRAS G12D | 2,632,299 | 258,897 | 258,897 | 254,364 | 4,533 | 3 |
| KRAS G12V | 4,080,670 | 1,979,594 | 1,000,000 | 983,854 | 16,146 | 10 |

These are per-target molecule IDs, not a cross-target unique union. G12V has a
quota-limited reserve. Source QC quarantines were not backfilled. Exported MOL2
records retain original coordinates and are not validated docking poses. The
exploratory grid profiles permit clashes and do not establish affinity/selectivity.

Evidence: `D:/agent/to_human/KRAS_E123_review_20261029/` (24 report/receipt
consistency checks) and `data/e124-core-integration/returned_report_check.json`
(new bridge import, expected counts, no queued computation). Remote MOL2 payload
verification and subsequent docking are outside these checks.

## Attach results in Chat

Keep the existing Chat storage root, runtime and model configuration. When a report
is outside the active Project, append explicit import directories to the usual
server command, for example on the school host:

```bash
--result-import-root "$HOME/kras-grid-production-v1" \
--result-import-root "$HOME/AF3_ligand"
```

These options allow bounded, read-only JSON report import. They do not authorize
arbitrary execution. Paths must exist on the Chat host; a local Chat process
cannot read school-host files by their Linux path. Copy a report to an allowed
local directory or run Chat on the host containing it.

In **Tasks & results > Attach an existing result**, choose:

- **KRAS pool export report**: final E123 `report.json`.
- **AF3 ligand batch report**: the batch's `run_.../report.json` or host-check report.
- **AF3 prediction manifest**: the prediction manifest created by the existing adapter.

Alternatively, use these deterministic commands without an LLM call. Substitute
the actual absolute path; JSON strings do not expand `$HOME` or `~`.

```text
/workflow_results {"operation":"attach","kind":"kras_export","report":"/home/y/yuzhang/kras-grid-production-v1/report.json"}
/workflow_results
/workflow_results {"operation":"check","attachment_id":"RETURNED_ID"}
/workflow_results {"operation":"read","attachment_id":"RETURNED_ID","pointer":"/exports/KRAS_G12D"}
```

Natural-language Chat also has the same read-only tool. Supply the exact path in
the request; the model must not invent one. An attachment is session-owned and
deduplicated by path/kind/content hash. Changing the source does not overwrite its
saved evidence: check it, then attach again to record a new version. A running
report is a point-in-time snapshot, not a live remote monitor.

## AF3 semantics and fixes

`af3_prepare.params.protein_copies` defaults to 1. Set 2 for an explicitly requested
homodimer. All copies use the same verified construct. Protein and ligand chain
IDs are allocated without collisions; literal stereochemistry, charge and seeds
are preserved. Preparation does not claim successful inference.

Multiple ligand entries in a single input describe one complex. For independent
candidate predictions, use separate prepare/run pairs or the standalone batch
runner. The optimized shared-MSA runner is not yet dispatched by Chat.

The batch runner now reports total jobs and records missing final models, invalid
confidence JSON and launch errors as failures rather than leaving a running
report or a misleading zero validation return code. The original process return
code is recorded separately. Host checks remain explicitly separate from inference.

## Local verification

Protocol: `experiments/E124-chat-core-integration.md`.

- AF3/plan suite: 31 tests passed after correcting the literal-SMILES test fixture.
- Result/domain/full-library/structure/budget suite: 56 tests passed.
- Final batch/result suite: 31 tests passed, including pending-job accounting.
- Existing Chat and batch tests also passed in the initial combined run.
- Actual returned KRAS report imported through Chat with expected counts and zero jobs.
- Frontend JavaScript syntax checked; authenticated HTTP attachment exercised.
- English-content guard passed for 1,006 maintained text files. The global diff
  check also reported pre-existing research-log whitespace; historical log lines
  were preserved rather than rewritten.

These suites overlap; their counts must not be summed as distinct tests. Final
checks and any limitations are recorded in the experiment log. No new library
scan, docking, GPU prediction or live-provider benchmark was run for this update.

## Remaining work, in order

1. Dispatch optimized independent-ligand AF3 batches through the Chat queue, with
   job-level progress, cancellation semantics and code/input-bound resume. Preserve
   fixed cofactors and distinguish independent ligands from co-ligands explicitly.
2. Add a typed Chat execution adapter for the latest retained-pool KRAS production
   plan, or a general equivalent with explicit target/query/receptor provenance.
   Do not present external result import as this execution adapter.
3. Validate both paths on the actual school runtime: container arguments, mounted
   databases, GPU execution, failed-job recovery and model/export payload checks.
4. Integrate imported prediction payloads with confidence review and PyMOL only
   after local paths and hashes are verified. Imported report evidence alone does
   not enable the existing `/confidence` or viewer task actions.
5. Revisit block routing after defining a matched endpoint and held-out utility
   evaluation. Keep it separate from reliable AF3/search infrastructure work.

Existing sealed execution plans remain tied to their code/configuration. This
update does not bypass resume mismatches or alter completed scientific artifacts.
