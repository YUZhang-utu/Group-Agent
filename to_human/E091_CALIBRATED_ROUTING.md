# Candidate-specific boundary placement with an explicit 1% target

The user confirmed a special-set target <=1%. The denominator is the admitted
24,663,736-conformer library, not the separate 1,150,072 chemistry-review records.
At most 246,637 admitted conformers can remain special to meet this target.
The <=45 angular candidates number 1,262,619; at least 1,137,702 must obtain
supported placement. This is a requirement to measure, not a reason to loosen
all thresholds. The preceding class-balanced panel cannot predict this count.

## What is implemented

- All <=35 and 35-45 candidates are evaluated individually, not only a pilot.
- Each candidate target gets up to 64 distinct-molecule references selected by
  a stable molecule-ID hash, rather than the first database entries.
- References split 50/25/25% into fit, calibration and independent check sets.
  A target needs at least 32 distinct references. Up to 12 diverse fit
  prototypes are selected, without using calibration or check molecules.
- Each comparison combines continuous side-chain properties, actual-coordinate
  steric proxies, and circular backbone torsion features. The routing backbone
  metric is sin/cos RMS, not the atom RMSD reported in E089.
- A joint score requires one prototype to support all three channels. It
  cannot use one prototype for chemistry and another for sterics to manufacture
  compatibility. The fixed calibration quantile is 95%; the independent-check
  gate is 80% retention. These are exploratory defaults, not activity thresholds.
- The chosen target, directed rotation, supporting prototype, normalized score,
  all accepted alternatives and original class are recorded for each CID.
- Sparse or failed models and property/steric/backbone failures remain special,
  with explicit reasons. No capacity-based splitting of logical classes occurs.

Reference check data are used to gate class models. Their measured retention
is therefore not a final independent estimate of whole-search recall after
selection. Future held-out retrieval and docking validation remain necessary.

## Run on the workstation

Update feature/structure-guided-chat. Run from the repository root in the
environment that completed the backbone build:

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BUILD=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone
python -m aidd_agent.calibrated_boundary_routing \
  --build "$BUILD" \
  --candidates "${BUILD}-boundary-candidates-v1" \
  --output "${BUILD}-calibrated-routing-v1" \
  --reference-limit 64
```

This is heavier than the 300-query panel: approximately 1.26M candidate
conformers plus references need source-verified chemistry/steric extraction.
It does not re-describe all 24M conformers, but source text scanning can be large.
It stores double-precision vectors in features.sqlite and logs ongoing progress.
No whole chemistry/typed build is required. Existing databases remain unchanged.
Use a fresh output; automatic interrupted-run resume is not implemented.

## Inspect the result

- report.json: actual regular/special counts, per-band counts, reason counts,
  whether <=1% is met, additional placements needed, and phase timings.
- boundary_assignments.sqlite: one row for every original boundary CID,
  assigned_provisional or special. Ordinary definite members are unchanged.
- routing_models.json: prototypes, calibration scales/cutoff, disjoint IDs and
  independent-check scores for every candidate target.
- reference_manifest.json: selected reference identities.
- features.sqlite: compact source-hashed properties/steric/backbone features.

The new assignment is a sidecar over the frozen source. It does not rewrite
old leaves or their stored radii. Old leaf bounds must NOT be reused as bounds
for enlarged logical classes. Production dispatch/bounds and incremental library
updates require a separate integration and validation step. No automatic block
rejection or claim of guaranteed recall is made.

If the target is missed, use reason counts: sparse references require additional
evidence or compatible alternative targets; steric failures require different
targets or a special subgroup, not indiscriminate angular relaxation. Report both
percentage and absolute special-set size. Special-set search timing is explicitly
not_run here; phase timings measure this build only.

## Validation

Eleven focused tests pass, including disjoint molecule sets, diverse prototypes,
joint-prototype support, directed rotations and failed-check rejection. Real local
small-file run conserves all 24 conformers, leaves the 21 definite members intact,
and keeps three boundary records special because reference sets are insufficient.
This deliberately does not force the 1% target on inadequate evidence. Input SQLite
hashes were unchanged. Whole-library counts remain pending workstation execution.
