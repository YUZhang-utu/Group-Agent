# Independent block search: restored consensus funnel

New Chat start/compare requests use `search_engine=consensus-threshold-v1`.
Existing requests retain `legacy-budget-v1`; they are not silently converted.

## What changed

- Block members enter the same threshold-mode compute function as the full-library consensus funnel. The adopted query supplies coarse constraints, anchor threshold, mandatory/alternative groups, geometry, pocket and final pose threshold. No new threshold or hard anchor is invented.
- Only actual per-template passes enter conformer ranking. Their union is ranked with the existing per-template contact/Gaussian ordering and conformer RRF. Whole-block boundary stopping remains: finish a started block, then stop if 100000 surviving conformers have accumulated. Shortfalls are not filled from rejected/outside-block candidates.
- Deterministic Hungarian column updates use NumPy, preserving reference traversal and first-tie assignment. There is no SciPy replacement solver, lowered seed cap or reduced template count.
- Completed block refinement is reused between independent arms only when query, code, input identities, exact member IDs, chunk schedule and backend match. Each arm still ranks its own searched prefix. Cached files are verified and local copies receive local hashes. Reuse is reported; its elapsed time is not an independent speed benchmark.
- Reports include engine, template IDs/count, filter policy and reuse timing scope. Progress includes current block size and template count.

The current saved MDM2 query has eleven co-crystal ligand templates. The previous handoff lists:
6Q9L:HTZ:A:201, 5C5A:NUT:A:201, 5LN2:6ZT:A:201, 7BJ6:TVK:A:207,
4OBA:2TW:A:501, 6Q96:HRE:A:201, 7BMG:U3Z:A:204, 7NA2:1I0:A:201,
3W69:LTZ:A:201, 3JZK:YIN:A:1, 5HMK:62Q:A:1001.
An earlier eight-ligand version has not been located. Do not remove three templates based on memory; a revised panel requires an explicit new adopted query.

## Safe workstation transition

Do not pull/restart while an old job is executing: protocols hash source code, and updated files must not mix with existing receipts. This update has not stopped task 7fa0b402ba4b4a3a or its queued peers.

When choosing to switch, cancel queued legacy tasks first using `/cancel TASK_ID`, then cancel the running task or let it complete. Confirm no active legacy computation remains, then stop the Chat server. Existing output files are retained. Do not resume legacy partial outputs with new code.

With the server stopped:

```bash
conda activate aidd-workstation
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git status --short
git pull --ff-only origin feature/structure-guided-chat
python scripts/check_block_assignment.py
```

If the checkout has local edits, preserve them and resolve the pull normally; do not reset them. Use the existing Chat launcher/runtime after the check. In the same conversation reuse the already imported EquiScore analysis and adopted query task IDs, and request a new search-only comparison. The new engine field creates distinct request identities; do not use `/resume` on the old tasks to migrate them. No scoring or docking should be requested.

Verify `/block_campaign {"operation":"comparison"}` reports the new engine. For a completed arm inspect `filter_policy`, `template_count`, `searched_conformers`, `ranked_conformers`, `shortfall`, `cache_reused` and `timing_scope`. More than 100000 input conformers may be needed, and a Top-5 arm may legitimately fall short.

## Evidence limits

Local fixture validation checks exact assignment indices (including ties and zero weights), inherited thresholds, full-funnel/explicit-ID equivalence, rejection before seeds, final threshold rejection, cache integrity and independent-arm stopping. The microbenchmark measures only synthetic assignment workloads. Actual workstation acceleration and current-target active retention remain to be measured. Optional-only designs may still admit many candidates; no rejection percentage is promised.

## User-authorized reset of the old twelve tasks

The user explicitly requested deleting the old comparison and starting again. Stop the existing Chat server with Ctrl+C first, and wait for its workers to exit. Then update the repository and run:

```bash
python scripts/clear_legacy_block_comparison.py \
  --storage-root /mnt/local/hand/yuzhang/aidd/e097-chat-workspace \
  --anchor-task 7fa0b402ba4b4a3a --apply
```

Omit --apply for an optional preview. The tool requires exactly twelve matching legacy search-only arms in the same session, validates sealed requests and bounded run paths, and refuses while Chat/compute processes remain. It removes those job rows, request files and run directories only. It retains a small removal inventory, conversation messages, upstream analysis/query files and other campaigns. It refuses ambiguous cohorts rather than guessing.

Restart with the established launcher, then ask Chat to compare E094/E095/E096 independently for ChemPLP and EquiScore, Top-5 and Top-10, using the existing adopted query and full analysis, stopping at 100000 surviving conformers and performing no docking. No new model scoring is needed.
