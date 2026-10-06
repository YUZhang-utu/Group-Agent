# AIDD-Agent continuation handoff - 2026-10-06

This operational handoff supersedes the 20261005 handoff for current progress.
Publication discussions are out of scope and must not be persisted as research
plans. Preserve the existing workstation task; inspect status before any submission.

## Latest status and immediate next action

The user reports the previously submitted leading-block search/docking task is
still running. Its NEW task ID and run path have not been supplied in this chat.
The exact current phase (mapping, search, export or docking) is not verified.
Do not claim 100000 candidates have already been produced or docking completed.

1. In the existing MDM2 conversation, inspect the task list or use
   `/block_campaign {"operation":"list"}` to locate the latest search task.
2. Read `/status NEW_TASK_ID`, `/results NEW_TASK_ID`, the task execution report,
   and `execution/campaign/blocks/progress.json` when present.
3. If running, monitor the existing task. Do not submit another search, restart
   its server, pull new scientific code, or change runtime inputs while it runs.
4. If complete, inspect actual searched/ranked/exported conformers, unique
   molecules, shortfall, scored pairs, failed jobs and output hashes. Review
   `docking_candidates.csv` and saved poses. Do not substitute planned counts.
5. If failed or blocked, inspect the precise error, stage outputs and receipts
   before choosing recovery. Do not delete sealed outputs or blindly resubmit.

There is no configured direct workstation connection in this assistant session.
Workstation execution status comes from user receipts, not local Windows tests.

## Intended submitted workload

- Partition/receptor: E095 / 7NA2_A.
- Leading blocks: union of EquiScore Top-5 and ChemPLP Top-5, deduplicated.
- Each scoring arm ranks blocks by the mean of its best 10 DISTINCT molecules,
  each represented by its best in-block scored conformer; not whole-block means.
- Reuse the confirmed 11-template MDM2 co-crystal LIGAND query.
- Search every conformer in the selected blocks, not the entire library.
- Retain at most 100000 conformers TOTAL across the union; multiple conformers
  of the same molecule are explicitly allowed. This is not 100000 unique molecules.
- Dock retained candidates against the existing reviewed 7NA2_A receptor/site.
- Report shortfall if fewer candidates are available. Do not pad outside blocks.
- A new EquiScore pass on these new poses is not automatically included.

The actual submitted request/selection remains the authoritative scope; inspect
it on the workstation because the new task receipt has not been returned here.

## Completed inputs (user workstation receipts)

Original project: `prj-adf8a1b9f4f8-prompt-aidd`.
Original execution root:

```text
/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution
```

- PLANTS: `plants-parallel20/report.json`, 8445 jobs, zero failures,
  821619 scored conformer/receptor pairs (273873 unique conformers, 3 receptors).
  Source hash: `e5a9487d23469ec4846021014a617ee5791af1f6688ec677b0509fe4ffd4dd39`.
- EquiScore full: `equiscore-v3/equiscore-full/report.json`, complete,
  821619/821619 pairs scored, zero failures, elapsed 50634.678872 seconds.
- EquiScore analysis: `equiscore-v3/equiscore-analysis/report.json`, complete,
  8283 ranked block/receptor entries, Top-10 distinct-molecule ranking.
- Successfully adopted into Chat, no new scoring:
  `PROMPT-2d0b02861b884ed5/execution/campaign/blocks/` under the same project's runs.
  Saved ranking SQLite/CSVs, candidates, comparison and review files exist per
  supplied receipt. Adoption full_pairs=821619, failed_pairs=0,
  computation_launched=false. The imported task ID was not pasted.
- ChemPLP is lower-better; EquiScore is higher-better. Do not average their units.
  Neither is measured affinity. Ligand chemistry review remains separate.

The saved query was reported complete, `kind=consensus_design`,
`readiness=ready_for_consensus_funnel`, with these ligand instances:

```text
6Q9L:HTZ:A:201   5C5A:NUT:A:201   5LN2:6ZT:A:201
7BJ6:TVK:A:207   4OBA:2TW:A:501   6Q96:HRE:A:201
7BMG:U3Z:A:204   7NA2:1I0:A:201   3W69:LTZ:A:201
3JZK:YIN:A:1     5HMK:62Q:A:1001
```

Discover its actual owned query task ID; do not invent one. Do not rebuild the
reference panel or invoke `/guided` for this block-restricted campaign.

## Configuration-block recovery already delivered

Old BLOCKED task: `f0656f4d2c584dc4`, not the new running task.
Old run: `PROMPT-f5ae731001a04619`.
Error: `Configure search.batch, block_evaluation.sampling_profile and plants_profile`.

Actual runtime file:

```text
/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/e097-config/runtime.json
```

- search.batch: `/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`
- search.workers: 16; coarse_chunk: 500; refine_chunk: 64;
  bounded_pair_seeds: false.
- sampling_profile and receptor_tools_profile: `sampling.json` and
  `receptor-tools.json` in the same e097-config directory.
- Missing separate plants_profile was the cause. The fixed implementation accepts
  the existing receptor_tools_profile for trusted PLANTS executable/scheduling,
  preserving original sealed receptor/site/search settings and engine hash checks.

Commits pushed:

- `f8bffb3`: recover original CLI co-crystal ligand queries in Chat.
- `e3dc314`: trusted tools-profile fallback and explicit `retry_config`.
- `6c11d12`: allow the known empty `campaign/blocks` executor directory during
  configuration retry; still reject scientific files or completed stage receipts.

The user was given this recovery command and later reported the task still running:

```text
/block_campaign {"operation":"retry_config","task_id":"f0656f4d2c584dc4"}
```

It creates one fresh sealed plan and preserves the blocked receipt. Repeating
the same retry reuses that fresh task. Do not issue it again merely to check status.
Do not use `/resume` across code/runtime changes. Local validation: 66 related
tests for the initial fix, then 39 campaign/workflow tests for the empty-directory
correction; English guard passed. These are software tests, not real-library proof.

## Environment and preservation

- Windows repository: `D:/agent/projects/aidd_structure_guided_chat`.
- Branch: `feature/structure-guided-chat`.
- Remote: `https://github.com/YUZhang-utu/Group-Agent.git`.
- Linux repository: `/mnt/medchem_taltio/wrk/yu_agent/Group-Agent`.
- Working EquiScore profile:
  `/mnt/local/hand/yuzhang/aidd/tools/equiscore-dl4s/profile.json`.
- Working environment: Torch 2.7.0+cu128, DGL 2.5.0+cu121, RDKit 2025.9.2,
  NumPy 1.26.4, MDAnalysis 2.7.0, ProLIF 1.1.0, RTX 5090.
  Preserve it; do not reinstall the incompatible old CUDA stack.
- Preserve unrelated dirty E103/20261005 handoff files and untracked review bundles.
- User communicates in Chinese; maintained UI/code/docs are English.
- The user wants Chat-based control and reuse of existing results, not manual
  rerunning of completed inference. Do not touch MOLIQ during this continuation.

## Resume prompt

Read AGENTS.md, research-state.yaml and
to_human/20261006_AIDD_AGENT_HANDOFF.md. Continue AIDD-Agent by locating the
already-running leading-block campaign in the original MDM2 conversation and
checking its current phase/receipts. Preserve its 100000-CONFORMER budget and
existing inputs. Do not resubmit search, old PLANTS or EquiScore jobs.
