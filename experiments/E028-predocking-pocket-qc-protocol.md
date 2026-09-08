# E028 Pre-docking aligned-pose and pocket-exclusion QC protocol

Status: protocol locked before implementation

## Fixed inputs

- Accepted E026 multi-cocrystal aggregation for 8BJU/QT9 and 1X8B/824.
- Accepted immutable conformer artifact catalog containing 299,999 conformers.
- Exact deposited receptor coordinate files corresponding to
  `8BJU-prepared-v1` and `1X8B-prepared-v1`.
- Primary pose: each docking task's retained
  `atomcentered_anchored_joint` candidate-to-query transform.

## Hypothesis

Applying each stored transform to the artifact-backed candidate heavy-atom
coordinates will reproduce a query-frame pose suitable for pre-docking pocket
QC. A coordinate-only protein exclusion check will identify obvious geometric
collisions before docking without changing retrieval membership or claiming a
chemically complete structure.

## Locked workflow

1. Verify aggregation outputs and conformer catalog inputs against hashes.
2. Select a deterministic Top-N task set independently for every query.
3. Resolve each stable global conformer ID from the immutable artifact catalog.
4. Apply the stored candidate-to-query 4x4 transform without optimization.
5. Read the exact query ligand and protein heavy atoms from the corresponding
   deposited mmCIF coordinate frame.
6. Measure candidate/query centroid distance, candidate-to-query nearest-point
   coverage, minimum candidate/protein distance, and counts/fractions below
   locked coordinate-only clash distances.
7. Export generic-element PDB point clouds for PyMOL inspection and JSONL QC
   records with source hashes and explicit representation limitations.
8. Emit a deterministic manifest and a PyMOL review script. Do not remove or
   reorder any retrieval or docking admission based on these exploratory QC
   values.

## Locked defaults

- Top tasks per query: 100.
- Query-neighborhood distance: 4.0 angstrom.
- Close protein distance: 2.0 angstrom.
- Severe protein distance: 1.5 angstrom.
- Candidate coordinates: artifact v1 heavy-atom centers after the stored
  transform.
- Protein coordinates: non-hydrogen polymer atoms; query ligand is excluded.

## Confirmatory acceptance criteria

- Stable global ID, molecule ID and conformer ID agree between every docking
  task and resolved artifact record.
- Every transform is finite, has shape 4x4 and has homogeneous last row
  `[0, 0, 0, 1]` within tolerance.
- Exported candidate coordinates equal direct transform application within
  floating-point tolerance.
- Every metric is finite and distance/count bounds are internally consistent.
- Exactly `min(Top-N, available tasks)` records are selected per query with
  deterministic tie handling.
- A repeated identical invocation produces identical JSONL and manifest content
  hashes.
- Retrieval, admission and docking-task source files remain unchanged.

## Representation and decision boundary

Artifact schema v1 stores heavy-atom coordinates but not atomic numbers, bond
orders or topology. E028 point-cloud PDB files therefore use generic `X` atoms
and are for geometry/PyMOL review only; they are not valid docking ligands and
must not be converted to SDF by guessing chemistry. Chemically correct SDF
export, projected-color features and terminal-torsion refinement require a
versioned, once-built artifact v2 companion. Coordinate-only clash metrics are
diagnostic annotations, not retrieval filters and not evidence of binding.
