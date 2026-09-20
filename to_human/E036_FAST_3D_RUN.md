# E036: parallel chunks and full-library 3D search

The user-reported E035 million-conformer run passed: 1,000,000 conformers,
952,402 source-grouped molecules, 6,000,000 comparisons, zero numerical error and
assignment differences, and exact conformer/molecule ranks. Remote artifacts were
not independently inspected here. Full E031 calls took 1.4476 seconds for QT9
(4.63% of historical refine time) and 1.9092 seconds for 824 (2.73%).
E031 remains annotation-only; 39/41 exported poses still require human review.

## Priority: actual full-library search

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e036_fast_3d.sh
```

Default output:
`/mnt/local/hand/yuzhang/aidd/e036-fast-3d/run-w16-c500-r64-v1`.
Default E034 source:
`/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118`.

The runner searches the existing whole-library FAISS index for two calibrated WEE1
queries at nprobe128/budget10000, recalculates exact candidate descriptor distances,
and runs Gaussian coarse/refine plus batched E031. This validates the execution
chain against saved outputs; it is not an accepted arbitrary-target service.

Scheduling uses 16 workers, coarse chunks of 500 and refine chunks of 64, versus
E034's 2000/250. Scoring, per-objective Top5000, the 512-seed cap and budget are
unchanged. Coarse tasks increase from 5 to 20; refinement has about 99-107 tasks.
This may reduce underutilization/tail waiting, but workstation timing is required.
At most 2 x workers futures are in flight. Queries run sequentially and share one
index load.

Every Gaussian array (including transforms) and E031 score/assignment/rank must
match E034. Fresh retrieved candidates must also match. Integrity/index loading
and final equivalence checks are separate; no full-library exact-reference scan
or million-conformer stress test is run. `execution_seconds` includes candidate
comparison and serialization; it is not an SLA.

Return `report.md` and `report.json`. Compare total query and coarse/refine wall
time, not just seconds per chunk. Real scheduling speedup remains unmeasured here.

```bash
bash scripts/run_e036_fast_3d.sh --resume
```

Completed queries retain original timings. Partial chunk reuse is marked
`fresh_compute:false` and cannot establish fresh-query latency. Use a new
`E036_OUTPUT` for a full timing run.

## Parallel stress validation (no need to repeat passed million evidence)

The E035 wrapper supports `E035_WORKERS` (default8) and `E035_CHUNK_SIZE` (default512).
Spawned workers own independent mmap readers; deterministic assembly and hash-bound
receipts remain. Logs show progress/ETA every10chunks; `scale/progress.json` updates
every chunk. ETA excludes final global ranking checks. Cumulative worker compute
and elapsed wall time are reported separately.

Optional fixed100k serial/parallel comparison on an idle workstation:

```bash
E035_WORKERS=1 E035_SCALE=100000 \
E035_OUTPUT=/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/scheduler-serial-100k-v2 \
bash scripts/run_e035_review_and_scale.sh

E035_WORKERS=8 E035_SCALE=100000 \
E035_OUTPUT=/mnt/local/hand/yuzhang/aidd/e035-review-and-scale/scheduler-parallel-100k-v2 \
bash scripts/run_e035_review_and_scale.sh
```

Compare `scale.execution.scale_wall_seconds_this_invocation`, requiring no reused
chunks, identical sample hashes and passing equivalence checks. Execution order
can affect caching; a single pair is exploratory. Eight workers is a conservative
starting point, not an eightfold speedup guarantee. Memory, I/O and molecule
complexity matter. Preserve `run-1m-v1`; new code must use new protocol directories.

If finer scheduling helps little, profile Gaussian seeds, geometry and reads.
Reducing candidates or seed caps is a separate coverage experiment. Molecule-level
Top100 overlap of6/9does not show E031 superiority; discordant poses need review.
