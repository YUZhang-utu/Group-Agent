# AIDD checkpoint: two workstation jobs running, 2026-10-01

## Immediate instruction

The user reports both the top-million conformer selector and E094 work-block job
are running on the workstation. Save continuity and wait for their outputs.
Do not start duplicate jobs, rerun search, pull into the running checkout, replace
standalone scripts, change parameters, or stop/restart either process. There is
no direct workstation process access. No completion receipt for either job has
been supplied. This supersedes the September 30 handoff's selector-not-started
next steps. Historical context: 20260930_AIDD_EXPORT_AND_BLOCKS_HANDOFF.md.

## Job 1: MDM2 top one million CONFORMERS

The prior all-selected export is user-reported complete: 1,000,000 molecules,
3,854,146 original conformers. Selection array lengths were pasted; its export
completion report has not been directly inspected here. The desired next result
is 1,000,000 conformers, allowing several conformers from the same molecule.

Standalone implementation: scripts/select_top_conformers.py, published commit
a6a3205, GitHub YUZhang-utu/Group-Agent, feature/structure-guided-chat. Downloaded
to /mnt/local/hand/yuzhang/aidd/export-tools, with export_selected_conformers.py
beside it. No workstation git pull was needed. Seven local related tests passed.

Search source:

```text
/mnt/local/hand/yuzhang/aidd/e054-chat-20260923-104806/users/workstation/projects/prj-42c81b1322ab-prompt-aidd/runs/PROMPT-59addb9ffd6444c3/execution/guided/guided/search
```

Exported original conformers:

```text
/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-1m-all-conformers
```

Running selection output, log and PID (as supplied in the launch instructions):

```text
/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers
/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers.log
/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers.pid
```

The old ranking.sqlite stores molecule/template representatives only. The new
selector reads the sealed search/chunks/TT/*.poses.jsonl and NPZ schedules to
reconstruct conformer-level ranks: contact first, Gaussian tie-break, then
saved-policy positive-contact tier / quota tier / RRF fusion. Hash, identity and
full schedule checks precede output. Missing scores are not replaced with source
order or molecule ranks. This is a new ranking over previously scored ANN
candidates, not exhaustive full-library scoring or a binding affinity estimate.

Expected outputs: selected-conformers.csv, selection.sqlite, top-00001.mol2
through top-00100.mol2 (10,000 each), per-part CSV/receipts, report.json. Coordinates
remain the original library coordinates, not aligned or docked poses.

Next: inspect the log or final report. Completion must show status complete,
input_conformers 3854146 and exported_conformers 1000000. Also review scored and
unscored counts and distinct_molecules; the latter may be below one million.
No runtime or finish time has been measured here. Read MDMM2 policy details in
MDM2_TOP_MILLION_CONFORMERS.md. If interrupted, only after confirming no active
process, repeat identical arguments with --resume. Never launch a second writer
against the same output. Completed template ranks and output parts are reused.

Earlier user reported an unexpected-indent error in the Python environment
check. A one-line shell python -c correction was supplied. The latest user says
the job is now running; do not treat the old copy/paste error as current failure.

## Job 2: E094 work-block finalization and side-chain properties

User reconfirms the earlier job is still running. Implementation commit 79b9d4e;
entry scripts/run_e094_work_blocks.sh all. Guide: E094_FINAL_WORK_BLOCKS.md.

Input:
/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone

Output:
/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/final-blocks-e094

Scope: 24,663,736 admitted conformers, ordinary work-block minimum 5,000, no
maximum capacity, preserved chemical identities within small-class execution
pools. Full side-chain profile coverage is required; refinement retains the
5,000 minimum and allows at most eight children per parent without forced splits.
E091 special pool is 132,885 (0.538787%). Separate chemistry-review population
1,150,072 remains deferred. E091 is finished; do not confuse it with running E094.

After completion inspect:

- backbone-validation.json: structural_gate passed.
- properties-validation.json: structural_gate and property_gate passed.
- properties/report.json: blocks_below_minimum zero, conserved populations,
  actual block numbers/sizes and property coverage.

Passing establishes offline block integrity, not search recall or docking
performance. Subsequent evaluation is 3D search, PLANTS/MDM2 and the user's N-E
workflow; the binding-affinity model is undecided. Both running jobs can compete
for disk/CPU resources; no new workload or performance intervention is requested.

## Suggested next-session prompt

Read AGENTS.md, research-state.yaml and
to_human/20261001_AIDD_TWO_RUNNING_JOBS_HANDOFF.md. Continue by reviewing the
MDM2 top-million conformer selection and E094 block-job logs/reports. Do not
rerun or modify either running job.
