# E037: one-shot workstation validation

Run this suite when the workstation is available. Separate E036 and repeated
million-conformer E035 runs are unnecessary.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e037_workstation_suite.sh
```

Defaults use the existing E034 and library directories. Output:

```text
/mnt/local/hand/yuzhang/aidd/e037-workstation-suite/run-v1
```

## Fixed comparisons

| Scenario | Coarse / refine chunk | Pair-seed generation |
|---|---|---|
| old-reference | 2000 / 250 | Generate all, then retain the prefix |
| fine-reference | 500 / 64 | Reference implementation |
| fine-bounded | 500 / 64 | Stop at 512 unique seeds |
| old-bounded | 2000 / 250 | Stop at 512 unique seeds |

Each scenario uses16workers and both calibrated WEE1 queries. Two rounds, the
second reversed, produce16timed query executions. A separate fine-bounded profiling
run covers both queries' first coarse/refine chunks and is excluded from latency
comparisons. The first chunk is a diagnostic sample, not every candidate's hotspot
profile. All scenarios retain the same first512unique seeds in the same order.

Each searches the existing full-library FAISS index, then runs Gaussian and optimized
E031. Budget10000, Top-N5000, nprobe128 and scoring remain fixed. Candidates must
exactly match E034; all Gaussian arrays/transforms and E031 scores/assignments/ranks
must pass equivalence checks. A failed scenario is logged while others continue;
the overall suite exits nonzero if any scenario fails.

## Outputs and resume

Return root `report.md` and `report.json`. The summary reports per-query medians,
speedups relative to old-reference, and timing eligibility. Reverse ordering reduces
but does not eliminate load/cache effects; two rounds remain exploratory, not an SLA.

`progress.json` tracks completed scenarios. Each has a root `.log` file; inspect the
active log with `tail -f` from another terminal. Child directories contain individual
reports and artifacts. Profiling files are under
`profile-fine-bounded/{8bju,1x8b}/gaussian/{coarse,refine}/`:
`worker-first-chunk-profile.txt` and `worker-first-chunk.pstats`.

```bash
bash scripts/run_e037_workstation_suite.sh --resume
```

Completed queries retain original measurements. Partial chunk reuse is excluded
from speed comparisons. Use a fresh `E037_OUTPUT` for new complete timings.
Code/parameter/input changes invalidate the old protocol; do not mix new methods
into existing output. Preserve E034, E035 and prior E036 evidence.

Overrides: `E037_WORKERS` (16), `E037_REPEATS` (2, minimum2), `E037_OUTPUT`,
`E037_E034_SOURCE`, `AIDD_BATCH`. Avoid concurrent heavy jobs while measuring.
Startup integrity checks and index loading are separate in child reports; there
is no full-library exact-reference scan.

## Historical local validation

At implementation:205tests passed,2dependency skips, including multiprocessing,
resume/tamper handling, failure continuation, reverse ordering, exact seed prefixes,
full Gaussian arrays and separate profiling. Twelve fixed-seed synthetic cases with
three repeats each preserved exact prefixes. Generator speedups were1.27-9.69xwhen
more than512seeds were generated; below the cap they were about0.99-1.01x.
Raw data: `E037_LOCAL_SEED_BENCHMARK.json`.

These are generator microbenchmarks, not full-query speedups. Real-library timing
and hotspots require the workstation suite. E035's39/41exported discordant poses
still need human review; this suite does not establish pose quality or enrichment.
