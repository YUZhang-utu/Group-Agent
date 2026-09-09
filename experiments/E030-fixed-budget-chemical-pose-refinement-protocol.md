# E030 Fixed-budget chemical pose refinement integration protocol

Status: pivoted before batch integration; rigid micro-seed primitive retained
for downstream docking pose preparation only

## 2026-09-09 scope correction

This protocol must not be integrated into the general pre-docking 3D search.
Thirteen rigid seeds times the bounded two-torsion beam can require roughly 403
state evaluations per pose, over three million for the accepted 7,704-member
rigid-refinement set. That violates the fast-search objective and overfits a
WEE1 pose-QC observation. The implemented seed generator has no production
caller and therefore has introduced zero search latency. E031 supersedes the
batch-integration portion and reserves this primitive for a later docking
adapter operating on already admitted tasks.

## Motivation and human evidence

Human review confirmed that all six E028 representatives occupy the intended
co-crystal pocket, ruling out a gross query/receptor frame error. For 8BJU/QT9,
the macrocycle cores do not penetrate the receptor and conflicts are mainly
terminal. The 1X8B/824 LOW pose is also terminal-conflict dominated, while its
MEDIAN and HIGH macrocycle cores sit too close to the receptor. Numeric vdW
penalty did not perfectly reproduce the visual MEDIAN/HIGH preference, so it
cannot become a sole ranking truth.

## Hypothesis

A fixed-budget post-rigid refinement that retains the zero-change pose, explores
at most two bounded terminal torsions, and optionally evaluates small rigid-body
perturbations will reduce element-aware protein penetration for terminal-conflict
poses without losing the Gaussian/pharmacophore evidence that admitted them.
8BJU should benefit mainly from torsions; some 1X8B poses should require both
torsion and rigid micro-relaxation.

## Immutable boundary

- FAISS, pharmacophore, Gaussian, molecule aggregation, and protected per-query
  admission remain unchanged.
- Only a caller-capped final detailed set is chemically refined.
- The baseline transform and coordinates are always retained byte-for-byte.
- No pose or molecule is removed by vdW or direction metrics.
- Outputs are pose hypotheses and annotations, not docking scores or affinities.

## Locked initial budget

- Maximum candidate poses: explicit caller Top-N.
- Maximum stored terminal torsions evaluated: 2, smallest moving side first.
- Torsion offsets: -60, -30, 0, +30, +60 degrees.
- Torsion beam width: 8; the all-zero state is mandatory.
- Rigid micro-seeds: identity; translations +/-0.25 angstrom on each Cartesian
  query-frame axis; rotations +/-5 degrees about each Cartesian query-frame
  axis. At most 13 rigid seeds total.
- Candidate directions rotate with rigid geometry; translations do not alter
  vectors.

## Stored pose variants

Retain named variants instead of collapsing uncalibrated metrics into one score:

1. `baseline`: exact incoming pose and scores.
2. `min_vdw`: lowest soft exclusion penalty in the bounded search.
3. `max_directional`: greatest valid projected-direction overlap.
4. `balanced`: Pareto-supported candidate using an explicitly stored,
   versioned normalization; exploratory and unable to alter admission.

Each variant stores its transform or coordinates, torsion angles, Gaussian
shape/color values under that exact pose, directional overlap, vdW primitives,
and whether rigid movement was used.

## Confirmatory evaluation on real E028 Top-200

- Reproduce all 200 baseline coordinates, IDs, query/receptor assignments, and
  baseline vdW metrics within numerical tolerance.
- Report per query fraction improved, vdW reduction, Gaussian/directional
  changes, runtime, and torsion-only versus rigid-plus-torsion winners.
- Reinspect the six locked LOW/MEDIAN/HIGH representatives.
- Report rank correlations and Top-K overlap for every exploratory variant;
  do not replace the accepted retrieval order.

## Acceptance criteria

- Zero-change state is always present and exactly reproducible.
- No invalid coordinates, bond-length changes outside numerical tolerance, or
  nondeterministic repeated output.
- At least one nonbaseline pose improves vdW on controlled terminal-conflict
  tests; negative real results remain valid evidence.
- Runtime and memory scale with final Top-N and fixed budgets, never with the
  299,999-conformer library.
- New tests cover identity preservation, deterministic seed enumeration,
  vector rotation, torsion bounds, and variant provenance.

## Decision boundary

If core penetration persists, route the pose to a real docking/local-minimization
adapter. Do not expand this micro-search until it resembles an unvalidated
docking engine. Hard vdW filters or final weights require redocking and held-out
enrichment calibration.
