# E047: avoid unnecessary whole-ligand scoring

## Workstation finding

The user-reported E046 GPU pilot passed all array checks on 256 spread IDs, but
all rows survived invariant conditions and none passed the final rule. NumPy with
22 workers took 1.586 s median; CuPy took 7.143 s median (4.868 s second repeat).
Installing a GPU library did not solve the full-library cost. See
`E046_WORKSTATION_GPU_RESULT.json`. Do not start another full scan on this evidence.

## New compute funnel

1. Every conformer enters the unchanged invariant feature tests.
2. Generate the exact original bounded seed set for each invariant survivor.
3. Evaluate only selected typed/directional anchor geometry at those seeds.
4. If no seed could meet the rule, skip whole-ligand Gaussian computation.
5. Otherwise, reuse the decoded conformer and seeds, and score ALL original seeds.
6. Apply the unchanged E031 assignment and explicit final condition at the original
   Gaussian-winning pose. Retain all passing molecules, without Top-K/Top-N.

This is not an extra chemical requirement. For each selected anchor, the largest
compatible spatial-times-angular score is an optimistic upper bound on its score
under any one-to-one assignment. If even that upper bound fails, the pose cannot
pass. ALL is evaluated inside one pose; hits from different poses cannot combine.
A conformer is skipped only if every existing seed fails. Near-threshold numerical
uncertainty is retained. Seed-generation policy, 512 pair-seed cap, original seed
ordering, shape/color objective, thresholds and selected anchors are unchanged.

Crucially, we do not discard individual nonpassing seeds from surviving molecules.
One may still win the Gaussian objective; removing it could select a different
pose and falsely turn a final rejection into a pass. Original winner competition
is preserved. Rejected rows contain explicit unscored masks and placeholder values,
not purported Gaussian scores. Conditional diagnostics still count actual passes.

## What remains expensive

Seed construction remains query-dependent and linear in the conformers passing
invariant tests. No claim of sublinear arbitrary-query retrieval or minute-scale
whole-library completion is made. A persistent index of the same broad two-anchor
necessary conditions would retain the same candidates; building one now would
add offline cost without demonstrated pruning. A future index must first show
selectivity and conservative bounds for the actual rule. Query-independent local
frames and cross-conformer GPU scheduling are follow-up candidates if measured
seed generation dominates after Gaussian work is removed.

Repeated rules can reuse a completed funnel through selection previews. Changing
conditions requires new coverage when previously skipped poses were not evaluated.
Library expansion requires coverage of new shards; a previous complete report must
never be represented as covering new conformers.

## Workstation validation

Pull the new code and restart the existing chat launcher with `--allow-compute`.
Use the same conversation and saved selection preview. Enter:

```text
/benchmark
```

Or:

```text
Benchmark the latest saved selection rule with the new pose-feasibility funnel.
Include spread-out library conformers, saved positive examples and poses near the
threshold. Compare final membership and surviving poses against the unfiltered
reference. Report seed-generation time, pose-feasibility rejection, Gaussian work
and elapsed throughput. Do not launch the full-library run or docking.
```

The pilot includes 256 spread IDs plus up to 32 saved positives and 32 near-threshold
rows from the classified sidecar, deduplicated. Saved poses may not reproduce their
labels under current reference computation; the report records the actual reference
positive count. If none pass, positive-membership validation is explicitly missing.
The enriched panel is a correctness stress panel, not an unbiased prevalence sample.

Compare the scenario with `pose_feasibility: true` to the NumPy baseline. Inspect:

- `prefilter_passed`: invariant survivors; this can still be the whole sample.
- `pose_feasibility_rejected`: additional conformers excluded before Gaussian work.
- `gaussian_evaluated`: actual expensive-scoring count.
- `worker_seconds.seed_generation_seconds` and `.pose_feasibility_seconds`:
  summed worker time, not wall latency.
- `checks.false_rejections`: must be zero; all surviving arrays must match exactly.
- `reference_positive_count`: positive examples actually reproduced by reference.
- `scalability_gate`: sampled speedup, reduced Gaussian work, and positive coverage.
  This is a pilot gate, not full-library scalability acceptance.

The GPU comparison remains available for continuity; CPU is the current measured
winner. Do not switch hardware just because a GPU is installed. Do not reuse old
outputs after code/protocol changes. Benchmark does not auto-launch full screening.

To apply the selected configuration before restarting chat, for example:

```bash
export AIDD_POSE_BACKEND=numpy
export AIDD_SCREEN_WORKERS=22
export AIDD_POSE_FEASIBILITY=1
```

`/funnel` then uses the latest saved selection rule. Set feasibility to `0` only
for a deliberate baseline comparison. Existing launcher thread limits remain.
New progress distinguishes invariant survivors from Gaussian candidates and final
matches. Local tests establish the implementation contract; real-library reduction
and elapsed-time benefit must be measured by the new workstation pilot.
