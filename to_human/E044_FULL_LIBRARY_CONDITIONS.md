# Uncapped full-library feature conditions

This is the expensive unfiltered baseline. The user-observed QT9 window processed
about 86 conformers/second. For a rule-based necessary-condition funnel, use the
[E045 guide](E045_NECESSARY_CONDITION_FUNNEL.md). Do not launch this baseline as a
substitute for a fast conditional screen.

E043 exhaustive retrieval compares every descriptor but only refines Top-10,000
conformers and a Gaussian Top-N union. It cannot count every library molecule
matching a feature condition. E044 is a separate, substantially heavier batch job:
every catalog conformer receives rigid-pose computation and classified feature
scoring, without either membership budget. It reuses the crystal query definition,
not the old retrieved/refined candidate IDs.

## Chat execution

Update the repository, stop the old chat server and restart it with the same
provider configuration and compute flag:

```bash
git pull --ff-only origin main
bash scripts/run_chat_agent.sh --allow-compute
```

In the conversation containing the completed classification task:

> Evaluate every library conformer against the features in the completed classification task. Do not use descriptor Top-K, Gaussian Top-N or a molecule cap. Save all feature scores and report full-library molecule counts at the existing diagnostic thresholds. Do not dock.

The deterministic equivalent is `/full_count CLASSIFICATION_TASK_ID`. This uses
the chat task ID, not the `PROMPT-...` directory identifier. `/full_count` without
an ID uses the latest task, which must be the completed classification. Do not
submit the old `search_3d` prompt to obtain uncapped condition counts.

All classified features are computed. The display levels 0.25, 0.5 and 0.75 are
exploratory diagnostics; they do not change computation membership. Original
E031 scores and rankings are not modified. Subsequent selections can use any
explicit score in (0,1], because raw feature scores are retained.

## Results and subsequent selection

`report.json` and `interactions.html` report evaluated conformers, distinct
source-grouped molecules, full coverage and per-feature diagnostic counts.
Completion requires every catalog ID and the accepted molecule total. A partial
or failed job never establishes full-library counts. Per-query counts are not
cross-query molecule unions.

After completion, use the new full-count task as the selection source:

> From full-count task FULL_COUNT_TASK_ID, preview molecules matching all of EXACT_ANCHOR_ID_1 and EXACT_ANCHOR_ID_2 in the same stored pose, with minimum score 0.5 and no molecule cap. Do not export yet.

Replace IDs with actual report entries; choose your own threshold. `any` is also
supported. The preview streams chunks and uses SQLite to deduplicate molecules,
selecting the highest Gaussian-score qualifying conformer, with stable global-ID
ties. Two different conformers cannot jointly satisfy a same-pose `all` rule.
An explicit export cap only limits exported representatives, never the reported
matching counts. Representatives are stored as JSONL in global-ID order to avoid
reopening random score chunks during export.

Then `/export SELECTION_TASK_ID` writes SDF and IDs; the existing pocket/docking
preparation can consume the handoff. Cross-query conditional unions/intersections
are not added by this change. Feature aliases share a score column, so multiple
crystal contact labels for one feature are not independent requirements.

## Scheduling, persistence and optional timing pilot

Chat uses eight persistent processes and chunks of 2,048 conformers. Chunk size
controls scheduling and file count only: nothing is selected or discarded within
a chunk. Scores, assignments, transforms and IDs remain in chunk NPZs. Counts use
a disk-backed molecule index. The existing query is small; the whole library is
not materialized into one global score array.

Use `/resume TASK_ID` after an interrupted job with unchanged code and inputs.
Verified chunks are reused; missing chunks are computed; derived counts are rebuilt
to avoid double counting. Changed code/protocol requires a new task. Integrity
preflight checks shard hashes; final integrity checks and receipt checks are real
costs. Resume timings are not fresh-compute timings.

For an optional explicit pilot using the user's completed classification:

```bash
python -m aidd_agent.full_library_screen \
  --classification /mnt/local/hand/yuzhang/aidd/prompt-workspace/users/workstation/projects/prj-252fa94197d6-prompt-aidd/runs/PROMPT-7976a65eadf54c9e/execution/screening/screening/report.json \
  --output /mnt/local/hand/yuzhang/aidd/e044-full-library/run1 \
  --workers 8 --chunk-size 2048 --max-chunks 8
```

This pilot is explicitly `partial`. Repeat the same command without `--max-chunks`
to process the rest of both query libraries and reuse completed chunks. Do not
run the CLI and chat against the same output simultaneously. Initial shard chunks
are not a representative performance benchmark; timing/storage projections from
them need qualification. Full-library runtime and storage have not been measured.

## Meaning of a match

The current algorithm selects one anchored-Gaussian objective pose per conformer
from principal-axis and bounded pair-alignment seeds, then evaluates feature
conditions on that pose. This preserves the existing pose protocol and feature
scoring. It does not exhaust all orientations, torsions, protonation states or
chemical interactions, and can miss another pose that would satisfy the features.
Thus the result is an uncapped count under this defined pose/scoring protocol,
not a proof that all possible matching molecules or biological actives were found.

No PLIP or docking is required. Halogen-specific features remain unsupported.
Diagnostic aggregation currently supports up to 63 distinct indexed features per
query, including all features in the two supplied WEE1 classifications; a larger
query fails explicitly rather than dropping features. Local tests cover coverage,
resume, cross-chunk deduplication, same-pose selection and Gaussian pose parity.
Real full-library chemistry, throughput and resource use require workstation testing.
