# E023 Staged, parallel, resumable Gaussian reranking protocol

Status: real staged workstation performance confirmed; interruption/resume and ranking analysis pending

## Motivation and observed baseline

The real workstation E022 run completed successfully on 1,000 QT9 candidates:

- pair-seed mode (`max_pair_seeds=512`): 219.70 seconds wall, 220.11 seconds
  user CPU, 100% CPU, 55,256 KiB maximum RSS;
- PCA-only mode (`max_pair_seeds=0`): 16.59 seconds wall, 17.97 seconds
  user CPU, 108% CPU, 54,916 KiB maximum RSS;
- both retained 1,000/1,000 IDs in input order;
- the admitted retrieval union contains 299,999 conformers.

The near equality of CPU and wall time demonstrates single-core compute
saturation, not storage or memory pressure. Applying the expensive pair search
to every admitted conformer is therefore the wrong resource allocation.

## Hypothesis

A deterministic two-stage pipeline can preserve the complete coarse-score
record while concentrating pair-seed work on the union of each objective's
Top-N. Process-level parallelism plus invariant self-overlap caching should
reduce wall time substantially without changing scores, selected IDs, or poses.
Atomic chunk checkpoints should make a repeated identical command resume after
power loss without recomputing valid completed chunks.

## Locked workflow

1. Coarse stage scores every input conformer with centroid/PCA seeds only.
2. Stable descending ranks are computed independently for shape-only,
   atom-centered unweighted joint, and atom-centered anchored joint.
3. Refinement input is the union of Top-N conformers from all three objectives;
   this is expensive-compute allocation, not candidate rejection. The complete
   coarse result remains immutable and available.
4. Refine stage scores that union with PCA plus typed pair seeds.
5. Candidate chunks are contiguous slices of the locked input order.
6. Each completed chunk is an atomic NPZ plus a hash-bearing manifest. A chunk
   is reusable only if configuration hash, ID slice, result hash, and output IDs
   all validate.
7. The run manifest locks all input hashes and scientific/numerical parameters.
   Reusing an output directory with different inputs or parameters fails.
8. Worker processes open their own read-only artifact maps and query package.
   Parent-process merge restores locked input order regardless of completion
   order.
9. Rerunning an identical command skips valid chunks; missing, partial, or
   corrupt chunks are recomputed atomically.

## Confirmatory acceptance criteria

- Optimized invariant-cached scores reproduce the E022 reference primitives,
  objective values, and best transforms within numerical tolerance.
- One-worker and two-worker runs produce identical merged arrays and selection.
- Coarse output contains every input ID exactly once in original order.
- Refine selection equals the deterministic union of per-objective Top-N and
  contains at most `3 * N` conformers.
- Refine output contains every selected ID exactly once.
- A valid pre-existing chunk is reused without changing its hash or timestamp.
- A corrupt/missing chunk is recomputed while other valid chunks are retained.
- Configuration mismatch in an existing run directory fails explicitly.
- Atomic completion leaves no `.partial` files.
- CLI supports `coarse`, `refine`, and `all`, plus workers, chunk size, Top-N,
  progress reporting, and resume control.
- Full dependency-light test suite remains green.

## Workstation confirmation plan

Run the same locked 1,000 candidates with one and 16 workers, compare all
arrays, then run the complete 299,999-conformer PCA coarse stage. Interrupt and
restart one real run to confirm checkpoint reuse before launching pair
refinement. Runtime scaling and final Top-N overlap are empirical results and
will not be inferred from offline synthetic tests.

## Offline results (2026-09-07)

- Implemented rigid-invariant query/candidate self-overlap caching; each seed
  now computes only the three required cross-overlaps.
- Implemented deterministic PCA coarse scoring over every input ID and stable
  per-objective Top-N union for pair refinement.
- Implemented process workers with per-worker read-only artifact/query state,
  contiguous chunks, progress/ETA reporting, and parent-order merge.
- Implemented atomic NPZ/JSON chunk pairs, configuration locking, checksum/ID
  validation, corrupt-chunk repair, and identical-command resume.
- One- and two-worker synthetic runs produced identical coarse, selection, and
  refine arrays. A deliberately corrupted chunk was recomputed while two valid
  chunks retained their hashes and timestamps. No partial files remained.
- Complete dependency-light suite passed: 92 tests.
- Real 1-vs-16-worker scaling, interruption recovery, full coarse runtime, and
  refined ranking remain workstation confirmation tasks.

## Real workstation results (2026-09-07)

- The 16-worker PCA coarse stage scored all 299,999 conformers in 40.83 seconds
  inside the runner (42.58 seconds whole-command wall), about 7,045 end-to-end
  conformers/s, while retaining every ID in 300 validated chunks.
- Stable per-objective Top-5,000 union selected 7,704 conformers, only 51.4% of
  the 15,000 no-overlap maximum, demonstrating substantial objective overlap.
- Pair refinement completed all 7,704 conformers in 77.99 seconds inside the
  runner (79.67 seconds whole-command wall), about 96.7 conformers/s.
- Refine used 591.86 user CPU seconds over 79.67 wall seconds (743% CPU) and
  171,356 KiB peak RSS with no major page faults or swaps.
- Only eight refine chunks existed at chunk size 1,000, so at most eight of the
  requested 16 workers could run concurrently; the observed 743% is consistent
  with saturating that available chunk parallelism rather than a worker fault.
- From-scratch coarse plus refine compute was about two minutes. Score/rank
  analysis and a real interrupted-command resume check remain pending.
