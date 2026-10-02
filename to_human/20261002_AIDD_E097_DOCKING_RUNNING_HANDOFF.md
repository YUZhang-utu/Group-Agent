# AIDD E097 docking in progress - 2026-10-02

This checkpoint supersedes the earlier same-day E097 implementation handoff.
Resume repository: D:/agent/projects/aidd_structure_guided_chat.
Remote: https://github.com/YUZhang-utu/Group-Agent.git.
Branch: feature/structure-guided-chat. Latest published code: 12e0412.

## Immediate operating rule

The user confirms docking is RUNNING on the workstation. Do not restart, resample,
reprepare, repartition, pull code into its running checkout, change sealed profiles,
or launch a second docking scheduler. First inspect existing process arguments and
progress. Last explicit numerical progress was 145/8445 under the original runner.
A switch from 4 to 20 workers was supplied, but the user has not confirmed the switch,
its PID, completed-checkpoint migration, or which directory is now active. The latest
message says only that docking is ongoing. Do not assume the 20-worker run exists.
No remote process was inspected or controlled by this assistant.

## Frozen inputs and completed preparation

- E094: 384 regular blocks; E095: 1188; E096: 1189.
- 100 conformers per regular block; seed 20261002; 276100 sampling slots.
- Shared panel: 273873 unique conformer CIDs. Same molecule's different conformers
  remain separate. Deduplication shares results only for the same receptor/settings.
- 2815 exported ligand batches times 3 receptors = 8445 docking jobs.
- Expected full successful panel: 821619 conformer/receptor scores.
- Special pool 132885 and unsupported chemistry are excluded, unchanged.
- Rigid ligand, ChemPLP, speed1; current retained clusters default 1, not the 5 shown
  in the user's example. No scoring or ligand chemistry changes were requested.
- PLANTS: /home/y/yuzhang/PLANTS1.2/PLANTS1.2_64bit.
- SPORES: /home/y/yuzhang/pdb/1/SPORES_64bit, mode complete.

The earlier 6Y4Q command example was subsequently used as the actual starting PDB in
the completed workstation target/receptor assessment. Do not repeat the old statement
that no actual reference was selected. Target/receptor assessment task ea9380bc47cc45f9
completed; its run is PROMPT-b3fd8d8eec214a93. User explicitly adopted three structures:

| Receptor | Selection ID | Reference ligand | Center (angstrom) | Radius |
| --- | --- | --- | --- | --- |
| 1RV1:B | pocket-f74ca11ba1da | 1RV1:IMZ:B:110 | 48.1062926829, 12.048, 34.6398048780 | 14 |
| 7NA2:A | pocket-a31d65095ba4 | 7NA2:1I0:A:201 | 2.5452564103, -7.3996666667, 7.6101538462 | 13 |
| 5J7F:A | pocket-7946f10a0b62 | 5J7F:6GG:A:201 | 45.3713404255, 7.0469361702, 40.0384468085 | 16 |

Receptor preparation task a55310a3a93e4479 completed: three SPORES receptors,
native reference-ligand heavy-atom geometric centers, padding 5 angstrom.
7NA2 and 5J7F have reported nearby removed other-chain components, including other
ligands approximately 2.9 angstrom away. Structural review of those contacts remains
pending; preparation success alone does not establish biological pocket suitability.

## Exact workstation paths and identifiers

Chat storage: /mnt/local/hand/yuzhang/aidd/e097-chat-workspace.
Session: 3e151854d74c4c179b1829f312fd6ea7.
Task DB: /mnt/local/hand/yuzhang/aidd/e097-chat-workspace/chat/conversations.sqlite3.
Agent traces live under chat/agents, NOT directly under storage/agents.
Runtime: /mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/e097-config/runtime.json.
LLM configuration directory: /mnt/local/hand/yuzhang/aidd/config/llm.
Chat port: 8766. Preserve workspace/session when restarting only after running work is understood.

Sampling/preparation task: 8d1383bfe8354ed0, stage sample_prepare, complete.
Define RUN as:

```bash
export RUN=/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution
```

- Outer execution report: $RUN/report.json (not a direct PLANTS input).
- Actual preparation report: $RUN/prepare/blocks/report.json.
- Actual sampling report: $RUN/sample/blocks/report.json.
- Original failed pilot, retain: $RUN/plants-first-batch/report.json.
- Verified recovered pilot, original full continuation: $RUN/plants-recovered/report.json.
- Original full continuation log: $RUN/plants-recovered-launch.log.
- Proposed 20-worker migration destination: $RUN/plants-parallel20/report.json.
- Proposed 20-worker log: $RUN/plants-parallel20-launch.log.

Full continuation was launched via CLI, NOT automatically registered as a new Chat
job. Do not tell the LLM to launch the same panel again. Existing preparation task
being complete does not mean CLI docking is complete. Agent import/adoption of this
CLI run remains future integration work.

## Real-engine evidence and fixes

First job 1RV1_B-000000 processed 100 conformers in 117.15 seconds, zero skipped.
PLANTS wrote rankings and poses, but old parser rejected Unknown PLANTS ranking header.
Native format has TOTAL_SCORE as the first header field while each data row has an
extra initial pose-name column. correspondingNames.csv and MOL2 titles match aliases.

Commit 863d9ad fixed this layout with strict row-width, finite-score, alias/entry and
pose checks, plus sealed fresh-output recovery without engine execution. User supplied
successful recovery: status partial, attempted_jobs 1, failed_jobs 0, scored 100,
first_job_gate passed. This is real first-receptor I/O acceptance, not cross-docking,
affinity or recall validation. Old pilot directory remains unchanged.

Commit c2210ce earlier improved planner feedback: sampling options must not be attached
to receptor preparation, detailed invalid-answer feedback, elapsed-time guard on invalid
retries. Live LLM quality is not guaranteed by local regression tests.

Commit 12e0412 adds --workers 1..64 for CLI scheduling and recover --include-checkpoints.
25 related local tests passed, including completed receipts newer than the aggregate
report, unchanged prepared inputs, worker-count changes and no redundant engine calls.
Source locking rejects migration from an active scheduler. A stopped old-version run
can be migrated to a fresh directory, importing verified completed receipts; unfinished
batches without receipts may repeat. A code-hash change forbids direct resume into the
old directory. Never edit a frozen profile or signature to bypass this.

20 workers was recommended for the user's 24-core workstation. Estimated full runtime
around 14 hours is a linear extrapolation from one job, not a measured parallel runtime.
Current worker count and actual utilization remain unconfirmed.

## Next-session read-only checks

1. Inspect process arguments before doing anything else:
   `ps -eo pid,args | grep '[a]idd_agent.block_plants run'`.
2. Identify the actual --output and --workers from the running command. If --workers
   is absent, read the sealed prepare/blocks/profile.json scheduling value (original 4).
3. Tail the matching launch log and read progress.json if present. A report.json can
   still be the old 100-score snapshot while the runner writes additional receipts.
4. Request only short fields, not enormous output_hashes: status, jobs, attempted_jobs,
   failed_jobs, scored, first_job_gate, workers_this_invocation. Do not claim completion
   until final report and process/log agree.
5. If complete, expected 8445 jobs / 821619 scores / 0 failures. Inspect failed job
   receipts if partial; never restart all work merely because the final report is partial.
6. Run block analysis against the actual finished run only after verification:

```bash
python -m aidd_agent.block_plants analyze \
  --source "$RUN/ACTUAL_DOCKING_DIRECTORY/report.json" \
  --output "$RUN/block-analysis-v1"
```

The placeholder must be replaced with the verified active/final directory. Analysis
outputs block_summary.csv, block_scores.csv and review.md, separated by scheme and
receptor. Missing results are not zero scores; population-weighted means are withheld
for incomplete panels. Preserve CID/receptor/pose-file/index for subsequent N-E rescoring.
No N-E or binding-affinity model has run. No docking-driven block rejection, search
recall, atomwise descriptor completeness, or superiority claim is justified yet.

The user asked only to save continuity this turn. No workstation actions were launched.
