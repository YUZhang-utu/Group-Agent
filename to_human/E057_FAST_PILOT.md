# E057 fast bounded consensus pilot

This is an exploratory throughput/selectivity probe of the current adopted contact
design, not an active/decoy validation or a promise of 100000 molecules in 30 minutes.
The command creates an adopted-design fixture with crystal controls, samples unique
registry molecule IDs uniformly using reservoir sampling, includes all their catalog
conformers, runs every selected template through production preselection compute,
persists real poses, unions molecules and computes the same scaffold keys.

No original library or recommendation is edited. This does not start the full scan.
Output must be new and outside the library. CPU thread counts should be fixed to
avoid nested BLAS oversubscription. Use the same runtime/library as Chat.

## First quick probe

Stop any concurrent screening workload before interpreting performance. After pull:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
conda activate aidd-workstation
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python -m aidd_agent.consensus_pilot \
  --recommendation /mnt/local/hand/yuzhang/aidd/protein-only-check-20260923-132230/proposal/report.json \
  --runtime /mnt/local/hand/yuzhang/aidd/config/prompt-runtime-apptainer.local.json \
  --output "/mnt/local/hand/yuzhang/aidd/pilot-256-$(date +%Y%m%d-%H%M%S)" \
  --molecules 256 --workers 24 --chunk-molecules 4 --seconds 600
```

The command uses the confirmed workstation runtime JSON. Alternatively use
`--batch /actual/library/batch` instead of `--runtime`. Never guess a different
library. The runtime reader accesses only search.batch and prints no credentials.

Sampling scans registry IDs and compact molecule-ID arrays, not all 3D coordinates.
This preparation is measured. The budget is checked during sampling, dispatch and
aggregation; an individual filesystem/SQLite call can overrun a check interval.
On scan timeout, owned workers are terminated and partial chunks retained. No
extrapolation is emitted for incomplete samples, because fast-finishing cases can
bias both output fractions and runtime estimates. Use a fresh output for another
probe; partial runs are not resumed as fresh timings.

## Interpret report.json

- `status`, `sample_complete`: only complete samples support extrapolation.
- `sampled_molecules`, `sampled_conformers`, `templates`: actual scope; molecules
  are library registry identities, not deduplicated scaffolds or conformers.
- `stages`: molecule survival after each stage, unioned across complete templates.
- `counts_summed_over_templates`: work counts, not unique molecule counts.
- `matching_molecules`, `pose_records`, `scaffold_groups`: distinct output volumes.
- `projected_matching_molecules` and `approximate_95_percent_matching_interval`:
  exploratory population estimate; zero sampled hits do not mean zero library hits.
  The Wilson interval does not include modeling error or guarantee biological recall.
- `preparation_seconds`, `scan_wall_seconds`, `merge_and_group_seconds` and
  `wall_seconds`: measured components including reads, worker initialization and
  writes. `summed_worker_seconds` identifies kernels; it is not elapsed wall time.
- `projected_processing_wall_hours` and `projected_worker_cpu_hours`: approximate
  scale-up under unchanged hardware/template/configuration. Full-library byte
  integrity, independent control retrieval and nonlinear merge/I/O effects are
  excluded. Repeated per-chunk initialization differs from the production scheduler.

Every chunk evaluates all templates; candidates are not truncated to fit output
capacity. Pose payloads retain template, mask and spatial signature. The pilot
does not implement the proposed dual-definition pocket scoring. It measures the
current protein-contact-only baseline; later scientific rule changes need a new run.

## Escalation to larger samples

Use 256 molecules for plumbing and gross latency only. If complete and fast, use
2000 in a new output, then consider 100000 with a 1800-second budget. Estimates
from tiny samples are especially uncertain for rare hits and heavy-tailed runtimes.
Report actual completion rather than promising that a chosen sample finishes.
Use a common saved proposal and RNG seed for comparisons. Do not choose thresholds
solely to achieve a desired output count; inspect selectivity and controls first.

Local tests cover reproducible molecule sampling, all-conformer inclusion across
shards, actual production compute on explicit IDs, multi-template orchestration,
timeout termination and suppression of incomplete-sample extrapolation. Synthetic
tests establish software behavior, not workstation throughput.

## Workstation resource policy

The workstation has 24 CPU cores and an RTX 5090. Start with 24 single-threaded
workers; do not multiply this by BLAS/OpenMP threads. Memory use and storage
contention can make fewer workers faster, so this is a measured starting point.
The current consensus preselection path uses CPU same-pose Gaussian scoring.
The separate full_library_screen CuPy backend does not accelerate this path;
setting AIDD_POSE_BACKEND=cupy would not enable GPU consensus screening.
GPU integration requires identical-score validation before production use.
The initial 256-molecule run has 64 chunks at four molecules per chunk.
It measures startup and scheduling as well as kernels, and is not a reliable
rare-hit population estimate. Local regression: 448 passed, 2 skipped.
