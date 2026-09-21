# E047: rule-aware rejection before expensive whole-shape scoring

The user-reported E046 pilot passed exact array comparisons on 256 spread IDs,
but all survived invariant bounds. NumPy/22: 161.38 rows/s; CuPy/1: 35.84 rows/s
including first use, 4.87 seconds on second repeat. Neither establishes acceptable
full-library scaling. Merely indexing the same nonselective bounds cannot help.

Hypothesis: testing selected typed, directional anchor geometry on the unchanged
ordered seed set is much cheaper than computing whole-shape Gaussian overlaps
for every seed. Reject a conformer only if NONE of its existing seeds could pass
the explicit rule. An optimistic independent-anchor maximum ignores competition
from other anchors, so any actual passing assignment must survive. For ALL,
conditions must be possible within the same seed; never OR across poses. Keep a
conservative numerical margin at the score threshold. On survival, score ALL
original seeds, not only plausible seeds: dropping a nonpassing but Gaussian-best
seed would change the original winner and final membership.

Reuse generated seeds and decoded artifacts for survivors. Record invariant
rejects, pose-feasibility rejects, seed work and Gaussian work separately. This
changes computed diagnostics for rejected conformers, not final rule membership.
Partial results are not full-library counts. No Top-K/Top-N is introduced.

Validation: deterministic no-hit, true-hit, threshold boundary, duplicate feature
assignment and cross-pose ALL cases; randomized comparisons against actual E031;
unchanged full scorer arrays using reused seeds. Benchmark both original and
optimized pipeline with identical membership and survivor poses/scores. Add saved
classified rows spanning positive/near-threshold evidence alongside spread IDs;
report missing positive evidence explicitly, never silently claim coverage.

This still has linear seed-generation work and is not a promise of sublinear
arbitrary-query search. Future persistent indexes must demonstrate selectivity
and lossless bounds for the actual rule before expensive library-wide builds.
Workstation pilot must measure rejection and wall time before another full scan.
