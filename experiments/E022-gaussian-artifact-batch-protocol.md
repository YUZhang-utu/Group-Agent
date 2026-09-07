# E022 Artifact-backed Gaussian batch reranking protocol

Status: offline confirmatory implementation passed; real QT9 workstation run pending

## Hypothesis

The immutable conformer artifact can feed the validated Gaussian reference
kernel directly, without reopening source MOL2 files, while preserving every L1
candidate. Deterministic invariant seeds should recover a known rigid placement,
and shape-only, unweighted joint, and anchor-weighted joint objectives should
retain independent optimal transforms.

## Scope and definitions

- Input identities are stable artifact `global_id` values.
- Coordinates are reconstructed as `origin + int16 / coordinate_scale`.
- Candidate heavy-atom and feature-center arrays are read from shard-local
  memory maps; source MOL2 files are not read.
- Query input is a versioned NPZ package containing heavy-atom coordinates,
  atom-centered feature coordinates/types, and anchor-derived feature weights.
- L1 admission is immutable: scoring may order candidates but never remove one.
- The three locked objectives are:
  - `shape_only`: shape query-biased Tversky;
  - `atomcentered_unweighted_joint`: mean of shape and atom-centered color
    query-biased Tversky with unit feature weights;
  - `atomcentered_anchored_joint`: the same joint score with query anchor
    weights.
- Each objective stores its own best transform and its complete primitive score
  record. Scores from different poses must not be mixed.
- Projected pharmacophore color is explicitly out of scope because artifact
  schema v1 stores feature centers but no directional projected points.

## Protocol

1. Resolve each global ID to exactly one contiguous artifact shard.
2. Validate shard schema, offsets, counts, ID consistency, and coordinate scale.
3. Reconstruct heavy atoms and feature centers from read-only memory maps.
4. Load and validate a portable query NPZ package plus its hash-bearing JSON
   manifest.
5. Generate compatible pharmacophore-pair seeds and deterministic principal-axis
   fallback seeds; include a centroid-translation identity seed.
6. Score every seed with the E021 kernel.
7. Select and retain the best transform independently for all three objectives.
8. Emit one compressed NPZ row per input ID plus a provenance/parameter manifest.

## Confirmatory acceptance criteria

- Synthetic artifact coordinates round-trip within the locked 0.01 angstrom
  quantization tolerance.
- Out-of-range, duplicate, fractional, or inconsistent global IDs fail clearly.
- Output IDs equal input IDs in the same order, including candidates with no
  compatible feature-pair seed.
- A known proper rigid transform is recovered by at least one deterministic seed
  on locked synthetic geometry.
- All three objective transforms are stored independently and their recorded
  scores reproduce when reevaluated at those transforms.
- Repeated runs with identical inputs are deterministic.
- Output manifests record formats, versions, operand direction, objective
  definitions, artifact/query hashes, parameters, counts, and limitations.
- The complete dependency-light test suite remains green.

## Later workstation experiment

Run the scorer on preregistered QT9 L1 IDs, inspect score distributions and
poses, and measure enrichment against the locked retrieval benchmark. No real
ranking, enrichment, or projected-color claim is made by the offline test.

## Offline results (2026-09-07)

- Implemented catalog-to-shard stable-ID resolution with relocation fallback,
  read-only mmap access, offset validation, and exact schema-v1 dequantization.
- Implemented hashed co-crystal query packaging and anchor-to-RDKit-feature
  mapping without reading the production library.
- Added centroid, principal-axis, and typed pharmacophore-pair seeds.
- Added batch scoring that retains every input ID and separately selects/stores
  shape-only, atom-centered unweighted joint, and atom-centered anchor-weighted
  joint poses with raw cross/self values.
- Locked synthetic geometry recovered the known rigid pose; stored scores
  reproduced at the stored transforms; reruns were numerically deterministic.
- Focused tests and the complete dependency-light suite passed: 90 tests.
- RDKit/Gemmi are absent from the laptop base interpreter, so real QT9 query
  packaging and artifact runtime/enrichment remain workstation experiments.
