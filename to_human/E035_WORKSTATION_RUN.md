# E035: molecule review and equivalent E031 acceleration

Reuse completed E034 outputs without rebuilding the library or rerunning Gaussian.
E031 remains annotation-only. The historical QT9 target was<=3.1242seconds for a
full call, not a kernel microbenchmark. Subsequent million-conformer workstation
results are in [the result record](E035_WORKSTATION_1M_RESULT.json).

## Run (default100000distinct conformers)

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e035_review_and_scale.sh
```

Default source:
`/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118`.
Current parallel-runner output:
`/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/run-100k-parallel-v2`.
The original serial handoff used `run-100k-v1`; preserve that historical directory.
See [E036](E036_FAST_3D_RUN.md) for8-worker scheduling and512-row chunks.

Use the existing Python/NumPy environment. Preserve full E034 evidence: report,
EXECUTION_COMPLETE marker, protocol, both original query/result NPZ files, manifests
and structures. A copied report alone is insufficient. Overrides:
`E035_E034_SOURCE`, `AIDD_BATCH`, `E035_OUTPUT`.

## Measurements

1. Verify report, query, rigid-result, sidecar and catalog provenance/hashes.
2. Profile the reference separately (`profile.txt`, `reference.pstats`), excluded
   from latency comparison.
3. For both queries, run three alternating reference/optimized full-call repeats,
   including reads, checks, all three objective poses, analysis and serialization.
   Compare every repeat to the original E034 sidecar.
4. Independently select each score's best conformer per source-grouped molecule;
   record IDs, cross-scores, ranks, correlations and Top-K overlap.
5. Select at most10molecules per group: Gaussian Top100/E031 rank>500, the reverse,
   and both Top100. Do not backfill empty groups. Export both winning conformers
   when distinct.
6. Sample100000distinct library conformers across shards. Two WEE1 anchor sets and
   three fixed centroid-aligned orientations give600000comparisons. Save sample IDs,
   shard coverage, errors and global rank checks. The original full readers verify
   every optimized-reader feature coordinate/type/direction/ID independently.

The optimized engine batches geometry and reads only required interaction features,
avoiding unrelated atom/bond/torsion objects. It retains the original Hungarian
assignment. Reference remains the default E031 engine; this experiment opts into
batched scoring for comparison.

Acceptance requires exact identities/order/assignments, absolute score and anchor
error<=1e-6, exact conformer ranks and molecule representatives/ranks. A tiny numerical
error that changes ordering fails. Failures stop and update `RUN_STATUS.json`.

## Pose review

Each query's `review/` contains `molecule-ranks.csv`, `poses.json` (IDs, transforms,
scores and anchor contributions), native-crystal-frame heavy-atom SDF files,
`receptor.cif`, and `review.py`. In PyMOL:

```text
run /absolute/path/to/review.py
```

Enable poses individually and inspect placement, receptor clashes, feature types,
directions and anchor contributions. Dashed assignment lines are not validated
hydrogen bonds. SDF exports do not reconstruct every original stereochemical label;
there is no new minimization or docking. Human review, redocking and enrichment are
separate validation tasks.

## Million-conformer tier

After the100kchecks, use a fresh output for the current code:

```bash
E035_SCALE=1000000 \
E035_OUTPUT=/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/run-1m-parallel-v2 \
bash scripts/run_e035_review_and_scale.sh
```

This is1Mdistinct conformers and6Mcomparisons, never duplicated old candidates.
The original serial `run-1m-v1` has already passed per user report; do not overwrite
or unnecessarily repeat it. Stress orientations are constructed, not1MGaussian
poses or biological validation. Original E034 pose checks cover6299+6791conformers
across three objectives, totaling39270scores.

## Resume and interpretation

```bash
bash scripts/run_e035_review_and_scale.sh --resume
```

Retain the same scale/output/worker/chunk settings. Inputs, code and environment
must match. Completed queries and512-row stress receipts are hash-checked. A partial
query timing stage reruns all repeats; stress resumes unfinished blocks. Existing
content without a protocol is rejected; absent output starts fresh.

Return root `report.md`/`report.json` and query `review/summary.json`. The gate uses
historical E034 refinement time with matching host/CPU count, not a newly measured
Gaussian denominator. Load/clocks are uncontrolled; three repeats are not an SLA.
Stress kernel timing excludes shared reads, checkpoints and global rank checks;
parallel cumulative worker time is not elapsed wall time.

Historical local evidence:100000synthetic kernel cases, zero assignment mismatches,
max error about1.1e-16, exactfloat32ranks, nonempty-group speedup2.1-3.7x;189tests
passed/2dependency skips. These did not establish workstation full-call performance.
The later real-library result is recorded separately rather than rewriting history.
