# Independent block-search comparison - 2026-10-07

The user confirmed whole-block sequential stopping, not exhaustive search over every Top-K block. Top-5/Top-10 denotes the maximum block count. The statistic ranking each block stays Top-10 DISTINCT molecules with one best in-block conformer per molecule.

## Workstation update

The prior union campaign completed according to the user receipt: run PROMPT-874365d177704bf3, E095/7NA2_A, 344476 searched/ranked conformers, 100000 exported/scored, 72859 distinct molecules, zero failed jobs. Preserve that run. Verify no other tasks are running before updating/restarting Chat. Do not pull new scientific code during an active run or resume old sealed runs across code changes.

From /mnt/medchem_taltio/wrk/yu_agent/Group-Agent, activate aidd-workstation and pull feature/structure-guided-chat. Restart Chat using the existing e097-chat-workspace storage root, runtime and provider configuration.

## Chat

In the original MDM2 conversation, inspect `/block_campaign {"operation":"list"}` to find the adopted full analysis task and confirmed query task. Natural language can request:

> Using the existing full EquiScore/ChemPLP analysis and confirmed MDM2 query, compare E094, E095 and E096 independently, with each scorer's Top-5 and Top-10 blocks ranked for 7NA2_A. Search one complete block at a time in rank order, stopping at 100000 conformers or the selected block limit. Do not dock. Keep MOL2 exports and return all task IDs.

Explicit command (replace both placeholders with actual owned task IDs):

```text
/block_campaign {"operation":"compare","task_id":"ANALYSIS_TASK_ID","query_task_id":"QUERY_TASK_ID","receptor":"7NA2_A","conformers":100000}
```

This queues twelve independent arms: three partitions x two scorers x two block limits. Repeating the same request reuses task IDs; failed/interrupted tasks are not blindly resubmitted. No union and no docking. Query/receptor ranking context is held fixed. If more than one confirmed query exists, select the original one explicitly.

Inspect all owned sequential arms:

```text
/block_campaign {"operation":"comparison"}
```

Inspect individual tasks with `/status TASK_ID` and `/results TASK_ID`. For a single arm, `start` accepts `method: "chemplp"` or `"equiscore"`, `blocks: 5` or `10`, `dock: false`, and `search_policy: "ranked_blocks_until_budget"` in addition to the analysis/query/scheme/receptor fields. Legacy start defaults remain unchanged; explicitly set dock=false for search-only work.

## Stopping and exports

Finish all conformers/templates of each started block. If cumulative ranked survivors reach the budget, stop before the next block, fuse the searched-prefix ranks with the existing contact-first conformer RRF, and export at most 100000. If Top-K runs out, retain the actual shortfall. Multiple conformers per molecule are allowed. Top-5 and Top-10 can be identical when both stop within the first five blocks. Rank truncation at the final block boundary is explicit; no claim to the optimum of unsearched blocks.

Each arm records searched_blocks, searched/ranked/exported_conformers, unique_molecules, shortfall, stop_reason, elapsed_seconds and resumed_execution. Member mapping may cover all selected blocks; searched_conformers counts only visited blocks. Resumed timing cannot establish fresh speedup. Arms currently run independently; shared-prefix compute is not cross-arm cached.

MOL2 files: execution/campaign/blocks/selected/ligands/*.mol2. They retain original source coordinates/types, with export aliases linked back by selected/exports.json and selected/candidates.jsonl. Search transforms remain in the candidate records. These are not docked poses. Source MOL2 chemistry is not repaired or inferred; protonation review remains separate.

Compare only matching query/analysis hashes and candidate budgets. Counts and retrieval ranks do not establish experimental enrichment or affinity. Current evidence is local fixture validation; the twelve real-library runs await workstation submission.
