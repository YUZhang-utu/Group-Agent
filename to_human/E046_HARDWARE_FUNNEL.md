# E046: diagnose the funnel and benchmark CPU / CUDA

## What the stopped run established

The user stopped E045 after the displayed 55,296 rows. Every displayed chunk
retained all 2,048 rows at the necessary-condition stage. Only a few final poses
passed. The prefilter was nonselective on that prefix; it did not prove that all
library molecules satisfy the rule. No completed full-library total exists.

Current search combines whole-ligand shape, typed atom-centered Gaussian features
and crystal-derived feature annotations. This is richer than a type-only feature
lookup, but superiority over a pharmacophore baseline has not been established.
The two selected hydrogen-bond hypotheses do not require hydrophobic, aromatic or
other classified hypotheses. Those only become requirements through an explicit
human rule. No additional anchor, threshold, ranking budget or docking contact
is silently imposed for performance.

## Changes

- Batch unchanged ordered rigid seeds in float64 and reuse color-distance kernels.
  Select potential winners before constructing full score records. Recompute all
  near-winning seeds with the reference CPU overlap kernel, preserving first-tie
  order. The 1e-10 objective contender guard is numerical protection, not a formal
  universal error bound for GPU arithmetic.
- Add conservative relative-direction pair bounds when directional query features
  are available. Vector and unoriented-axis kinds are handled separately. Broad
  ANY / two-feature rules can still retain most of the library.
- CPU defaults to affinity-visible logical CPUs minus two, capped at 24; 24 visible
  CPUs therefore select 22 workers. Full-screen chunks default to 256 rather than
  2,048 for earlier feedback and better load distribution. This scheduling choice
  is not a measured optimum for every library.
- Optional CuPy float64 overlap batches have one GPU owner process. Seed generation,
  reference winner checks and E031 remain on CPU. This is seed batching per
  conformer, not yet a cross-conformer GPU pipeline. A GPU may be slower; the
  benchmark decides empirically. No 5090 speedup has been measured locally.
- Receipts record prefilter, pose and annotation time. Progress records survival
  and throughput, and warns when early necessary-condition rejection is zero.

All library IDs still participate in a full funnel. Every survivor is precisely
scored under the existing heuristic rigid seed protocol. No Top-K or Top-N cap is
introduced. The protocol does not search every possible orientation or torsion.
Classification aliases for one ligand feature are not independent evidence.
E031 remains annotation-only in the original ranking; explicit user rules define
a separate membership selection.

## Workstation: small benchmark first

Update with `git pull origin main`, then restart the chat agent using the existing
launch command and `--allow-compute`. Existing sessions and saved selection rules
remain the source. Do not resume an old output after a code change; use a new task.

In the same conversation, enter:

```text
Benchmark the latest saved selection rule on a small sample spread across the
whole library. Compare the reference CPU implementation, batched CPU execution
and the GPU if available. Check scores, winning poses, E031 assignments and final
membership. Report prefilter survival, timings and the fastest equivalent tested
configuration. Do not start full-library screening or docking.
```

Deterministic equivalent:

```text
/benchmark
```

With no ID, this finds the latest completed selection preview in the current
conversation, even if the most recent task was the stopped funnel. The response
names the chosen source task. An explicit `/benchmark TASK_ID` accepts the saved
selection-preview task ID. This task never automatically starts the full run.

The default pilot uses 256 evenly spread global IDs, two repeats, a reference CPU
pool of up to eight workers, the same-sized batched CPU pool, the affinity-sized
batched CPU pool, and optional CuPy with one worker. Every sampled conformer is
scored, including prefilter rejects, to audit false rejection. Reports compare all
Gaussian arrays exactly, including winning seed IDs/transforms, and all E031
scores/assignments/condition decisions. A mismatch fails acceptance. GPU absence
is reported as unavailable, never as a passed GPU test.

Timing includes sample reads, process startup and result IPC. The first GPU repeat
can include kernel compilation; median of two is descriptive only. Worker-stage
seconds are summed CPU work, not wall latency. This pilot has no full raw-library
integrity sweep and no full output-serialization workload. It cannot establish a
service SLA, complete-library runtime, or biological recall.

For a larger independent pilot, use a fresh directory:

```bash
python -m aidd_agent.funnel_benchmark   --selection /absolute/path/to/selection/screening/report.json   --output /mnt/local/hand/yuzhang/aidd/e046-hardware/pilot-4096   --count 4096 --repeats 3
```

## GPU environment

The AF3 container's CUDA environment does not install CuPy in the host
`aidd-workstation` environment. Inspect the host environment before installing:

```bash
nvidia-smi
python -c "import cupy as cp; cp.show_config(); print(cp.arange(16).sum().item())"
```

If CuPy is absent, use the wheel appropriate to the host CUDA environment. For a
compatible CUDA 12 environment with runtime components supplied by wheels:

```bash
python -m pip install 'cupy-cuda12x[ctk]'
```

Use only one CuPy distribution. Do not install both cuda12x and cuda13x wheels.
The official installation guide documents the CUDA-specific packages and runtime
options: https://docs.cupy.dev/en/stable/install.html . Performance timing must
synchronize GPU work: https://docs.cupy.dev/en/stable/user_guide/performance.html .
This backend returns overlaps to CPU each batch, synchronizing the measured work.

## Apply a measured configuration

Set variables in the terminal BEFORE restarting the chat agent. Example for a
CPU configuration, with worker count taken from the benchmark recommendation:

```bash
export AIDD_POSE_BACKEND=numpy
export AIDD_SCREEN_WORKERS=22
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Only if CUDA passes equivalence and is faster on the workstation pilot:

```bash
export AIDD_POSE_BACKEND=cupy
export AIDD_SCREEN_WORKERS=1
```

Keep all existing launcher arguments and profiles. Then `/funnel` uses the latest
completed selection preview and the chosen local hardware settings. `/benchmark`
does not change these settings automatically. Backend is pinned in the new run
protocol. Source query/selection and integrity checks remain required.

## Local evidence

The local five-conformer synthetic CPU fixture produced exactly equal complete
Gaussian arrays and about 4.04x median scoring speedup over three repeats.
See `E046_LOCAL_CPU_BENCHMARK.json` for raw times and scope. It excludes library
I/O and E031, and is not a workstation or full-library speedup claim. CUDA tests
are skipped when CuPy is absent. Real pruning, CUDA equivalence and throughput
remain workstation validation tasks.
