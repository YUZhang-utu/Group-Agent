# AIDD / MEDCHEM Agent continuation record

Saved 2026-10-05 at the user's request. This is the current continuation entry;
older pending statuses elsewhere may describe earlier stages of the project.

## Current position

EquiScore is the user-selected public rescoring model. The real workstation
technical pilot completed: **96 pairs / 96 scored / 0 failed**, status complete,
scope pilot, exit code 0, elapsed 12.190821549855173 seconds. The user pasted the
report; the assistant has not directly read the remote score table. Local worker
and coordinator hashes match the pasted receipt. This proves technical execution
on this sample, not predictive accuracy, affinity calibration or full-panel success.

The assistant provided a background command for full scoring followed by analysis.
**The user has not confirmed launching it or supplied a full-run receipt.** First
inspect the workstation's existing process/log/report state; do not submit a
duplicate run. Preserve the successful pilot's code, profile and environment.

## Repository and accepted workflow

- Windows repository: `D:/agent/projects/aidd_structure_guided_chat`.
- Workstation repository: `/mnt/medchem_taltio/wrk/yu_agent/Group-Agent`.
- Remote: `https://github.com/YUZhang-utu/Group-Agent.git`.
- Branch: `feature/structure-guided-chat`.
- Latest executable-code commit: `9b917a1` (chemistry ABI correction).
- User accepted the MEDCHEM Agent Chat design, molecular background, Helvetica
  typography and hidden workflow/task panels in earlier conversation.
- Existing terminal PLANTS output has been adopted/analyzed through Chat; the
  previous unknown-tool/arguments problem was repaired. Use `/dock_results` to
  inspect existing attachment and task receipts before submitting anything.
- EquiScore terminal execution is implemented. **EquiScore natural-language Chat
  dispatch is not registered yet.** Do not route a rescoring request to docking.

Primary block statistic: mean of the **best 10 distinct molecules**, each using
its best in-block conformer, separately within scheme and receptor. Whole-block
means are diagnostic only. Top-5 block/candidate views are required. ChemPLP is
lower-is-better; EquiScore is higher-is-better. Reselect the top molecules from
all scored poses independently for each model, rather than rescoring only the
ChemPLP shortlist. Failed/missing values are never filled with zero. Incomplete
or undersized panels do not receive competitive primary ranks.

## Completed PLANTS and ChemPLP evidence

User-reported completed run: 8445 jobs, 8445 attempted, 0 failed, 821619 scored
conformer/receptor pairs, first-job gate passed. 273873 distinct conformer IDs,
three MDM2 receptors: `1RV1_B`, `7NA2_A`, `5J7F_A`.

Project: `prj-adf8a1b9f4f8-prompt-aidd`.
Run: `PROMPT-1d4f82d50d294925`.
Execution root:

```text
/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution
```

Relative paths: `plants-parallel20/report.json`, `plants-parallel20/scores.csv`,
`prepare/blocks/report.json`, `sample/blocks/report.json`.
Original docking report SHA256:
`e5a9487d23469ec4846021014a617ee5791af1f6688ec677b0509fe4ffd4dd39`.

ChemPLP Top-10 analysis completed per user receipt: 8283 scheme/block/receptor
records = (384 E094 + 1188 E095 + 1189 E096 blocks) times three receptors.
E095 / 7NA2_A leading block: `work-87f26357d1177ffd794707d5`, mean -93.17528.
User also retrieved its five distinct molecules and saved poses. Original ligand
mode is rigid; `ligand_chemistry_reviewed` remains false. No new docking needed.

## Working EquiScore environment

- Base interpreter: `/home/y/yuzhang/envs/dl4s/bin/python`.
- EquiScore worker interpreter:
  `/mnt/local/hand/yuzhang/aidd/tools/equiscore-dl4s/env/bin/python`.
- Profile: `/mnt/local/hand/yuzhang/aidd/tools/equiscore-dl4s/profile.json`.
- Official code: `/mnt/local/hand/yuzhang/aidd/tools/equiscore/EquiScore`.
- GPU: NVIDIA GeForce RTX 5090.
- Pilot: Python 3.9.25, Torch 2.7.0+cu128, DGL 2.5.0+cu121, RDKit 2025.9.2,
  NumPy 1.26.4, MDAnalysis 2.7.0, ProLIF 1.1.0.

The dedicated venv inherits base Torch/DGL/RDKit; its pinned chemistry packages
shadow base versions without modifying dl4s. The existing installer is now
`scripts/setup_equiscore_existing.py`. **Do not rerun installation for this
already-successful pilot or revert to the legacy Torch 1.11/CUDA 11.3 recipe.**
The coordinator still runs under the working AIDD Python 3.11+ launcher environment;
the profile selects Python 3.9 for the worker. Do not confuse the two interpreters.

Official commit: `8b2a9289cf7d181fa49de6ac6712260e8c500c4a`.
Checkpoint: `workdir/official_weight/save_model_screen.pt`.
Weight SHA256: `d4367bb73686b2363e238abb778fab55e2924458ec1ced561072bd82f711695d`.
Worker SHA256: `80b3627330752a3ae512d92deb206b953cc6bda01df4166bbc290d45a7bcdd47`.
Coordinator SHA256: `9b4ff7a5933d4b6bc336602bb6b6d82a59b8d3d67ff14bd28677d1c315d1ba92`.
Profile SHA256: `a04dd685614a28426e46d36d48af1aefc5d3fe3c20d15bc2e462a68d9c92927c`.
Pilot scores SHA256: `26945b9ce30e3bf70716c0073a464ed980d918f178ff014ea9a03dabfb89b71c`.

Fixed issues: known unused `mu`/`dev` checkpoint keys; hash-gated legacy checkpoint
loading under Torch 2.7; DGL historical function alias; retaining venv interpreter
symlink paths; NumPy/MDAnalysis ABI mismatch; active overlay package receipts.
Actual Torch/DGL kernels and an official-example full forward pass precede pilot.
Local final focused tests: 33 passed. Earlier full regression: 780 passed, 6 skipped
(before later compatibility patches); do not describe that as a rerun of final code.

## First actions next session

Set paths for read-only checks in the workstation repository:

```bash
export EQUISCORE_PROFILE=/mnt/local/hand/yuzhang/aidd/tools/equiscore-dl4s/profile.json
export EQUISCORE_OUTPUT=/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/equiscore-v3
pgrep -af 'aidd_agent.equiscore|equiscore_worker.py|run_e103_equiscore.sh'
tail -n 30 "$EQUISCORE_OUTPUT/full-and-analysis.log"
tail -n 30 "$EQUISCORE_OUTPUT/equiscore-full/worker.log"
```

Missing log files mean that those files are absent, not proof the job never ran.
Inspect `equiscore-full/report.json` and `equiscore-analysis/report.json` if present.
The success target for the full panel is complete / 821619 pairs / 821619 scored /
0 failed; validate actual receipts instead of assuming success from exit code.

The already-provided command (reference only, **not automatically resubmit**):

```bash
nohup bash -c 'bash scripts/run_e103_equiscore.sh full && bash scripts/run_e103_equiscore.sh analyze' \
  > "$EQUISCORE_OUTPUT/full-and-analysis.log" 2>&1 < /dev/null &
```

Interrupted inference supports same-input resume; explicit failed pairs remain
recorded and are not silently retried. Analysis requires a fresh output directory.
If full inference completed but analysis did not start, invoke only analysis after
checking for an existing process/output. Never overwrite a running/sealed analysis.

Analysis outputs under `equiscore-v3/equiscore-analysis/`:
`block_rankings.csv`, `block_top_candidates.csv`, `rank_comparison.csv`,
`block_ranking.sqlite`, `ranking_review.md`, `chemplp_baseline/`, `report.json`.
Compare ranks and top-molecule overlap per scheme/receptor. Inspect leading blocks'
Top-5 molecules and original pose files. Do not average ChemPLP and EquiScore units.

Then implement owned-result adoption and EquiScore dispatch/results in Chat,
respecting session/task ownership, input seals, pilot/full gates and idempotency.
Keep the N-E adapter interface for the user's forthcoming model code. Experimental
labels are not yet available/confirmed: EF1%, BEDROC and affinity validation remain
pending, and model agreement must not be labelled experimental enrichment.

## Resume instruction

Read AGENTS.md, research-state.yaml and this handoff first, then findings.md,
research_log.md, experiments/E103-equiscore-rescoring.md and
to_human/20261005_EQUISCORE_SETUP.md. Continue from the successful 96/96 pilot.
First establish the status of full scoring and analysis; preserve the installed
environment and existing outputs. Do not repeat docking or installation.
