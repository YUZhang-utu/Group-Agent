# Full-library necessary-condition funnel

The E044 unfiltered run is a baseline, not the intended fast funnel. User logs
showed 36,864 additional QT9 conformers in 429 seconds (85.9/s). Extrapolating
that short interval to the remaining QT9 IDs gives about 83 hours, excluding
other queries and aggregation. This is an uncertain local-rate estimate, not an
end-to-end measurement. Stop the old task with `/cancel TASK_ID` and keep its files.

## Membership contract

Every library conformer participates. Before expensive pose scoring:

1. Check required feature types and direction kinds.
2. Check that distinct query features can map to distinct candidate features.
3. For an ALL rule, check rigid-invariant pair-distance necessary bounds and
   repeatedly remove unsupported assignments.
4. Compute the unchanged Gaussian pose protocol and feature scores for every
   surviving conformer, then test the actual same-pose rule.
5. Retain all passing molecules, deduplicated across conformers. No Top-K or
   Gaussian Top-N limits membership.

For the existing score `exp(-d^2/(2*sigma^2)) * angular`, angular is at most one.
A required score at least t therefore requires d <= sigma*sqrt(-2*log(t)), also
bounded by the scoring cutoff. Pair-distance discrepancy cannot exceed the sum
of the two positional radii. The implementation expands the bound with a small
numerical margin, retaining uncertain boundary cases. Passing these tests is not
proof of a matching pose; it only means the expensive stage is still necessary.

This preserves the final pass set under the existing one-pose-per-conformer
protocol. That pose protocol itself is heuristic, not an enumeration of all
orientations or torsions. Biology and docking quality remain unvalidated.

## Use an explicit rule, not an implicit change of question

Start from the completed classification. Save the desired exact anchor IDs,
ALL/ANY and threshold as a selection preview using the existing chat interaction.
The preview's counts describe its source candidate set. The new funnel uses only
its rule and crystal query; it does not restrict the library to those candidate
IDs, and ignores any preview export cap. Even a zero-hit preview can be a valid
full-library rule.

Then send:

> Apply the rule from selection task SELECTION_TASK_ID to the entire library using the necessary-condition funnel. Reject only conformers that cannot satisfy the rule, refine all survivors, and retain all passing molecules without Top-K or Top-N. Do not dock.

Or `/funnel SELECTION_TASK_ID`. Update code and restart the chat server with
`--allow-compute` before creating the new task. Do not resume an old code-bound
E044 task with changed code. The full E044 run is not a prerequisite.

The rule must already be explicit: listing all interaction classes does not mean
every class is mandatory. No default threshold or required anchor combination is
invented. An ANY rule or separate single-feature diagnostic counts can usually
only use feature compatibility before pose computation, so substantial speedup
is not guaranteed. A selective ALL rule may admit stronger geometric rejection.
The implementation never tightens a rule just to make it faster.

## Outputs and subsequent selection

`progress.json` reports actual completed row counts, not the largest ID of an
out-of-order worker result. Logs show prefilter retention and matching counts.
The final report separates all checked conformers, prefilter rejects, survivors,
pose-evaluated conformers, matching conformers and distinct matching molecules.
Unscored rows have `pose_evaluated=false` and assignment value -2. Their zero scores
and identity transforms are placeholders, not calculated poses.

The HTML's per-feature counts concern only conformers passing the chosen rule.
They must not be presented as unconditional full-library feature totals. Reusing
a previously evaluated chunk can make pose-evaluated rows exceed prefilter
survivors; the receipt labels that reuse explicitly.

To export, preview the SAME rule from the completed funnel task, then separately
request `/export SELECTION_TASK_ID`. Other conditions require a new funnel run:
poses rejected by the old rule were not computed for arbitrary later conditions.
Cross-query conditional unions/intersections remain a separate capability.

## Optional pilot and reuse of the stopped E044 computation

The CLI can reuse compatible completed E044 blocks while applying and auditing
the prefilter. Supply the new selection-preview report path printed by chat:

```bash
python -m aidd_agent.full_library_screen \
  --selection /absolute/path/to/selection/screening/report.json \
  --output /mnt/local/hand/yuzhang/aidd/e045-condition-funnel/run1 \
  --reuse /mnt/local/hand/yuzhang/aidd/prompt-workspace/users/workstation/projects/prj-252fa94197d6-prompt-aidd/runs/PROMPT-1861eb9c89504b31/execution/screening/screening \
  --workers 8 --chunk-size 2048 --max-chunks 8
```

This is a partial pilot. It can measure necessary-condition survival and audit
existing passing poses, but reused blocks do not measure fresh pose-computation
speed. Remove `--max-chunks` to continue all IDs with unchanged parameters and
code. Reuse requires matching scientific code, NumPy, classification, query,
catalogs, pose parameters, chunk size and verified chunk receipts. A contradiction
between a passing legacy pose and the prefilter is a hard failure. A rejected
reuse attempt never silently changes the scoring protocol.

The prefilter has analytic necessary-condition justification and local boundary,
random-transform, assignment and pipeline tests. Real-library survivor rates,
overall speedup and resource use remain unmeasured for this mode. Do not promise
seconds or minutes until those measurements exist.
