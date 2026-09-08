# E029 One-time chemical/directional companion and fast query refinement protocol

Status: protocol locked before implementation

## Fixed architecture

- Existing artifact v1 coordinates, pharmacophore postings and FAISS indices
  remain immutable and continue to define protected broad recall.
- A versioned companion is built once per library shard during one sequential
  source-MOL2 pass. It is append-only and keyed by the same stable global IDs.
- Directional color, protein exclusion and limited torsion refinement execute
  only after fixed-budget rigid Gaussian refinement; they never scan the full
  library and never remove stored retrieval evidence.

## Hypothesis

Persisting chemical topology and conformer-local directional/torsion metadata
once will make chemically aware post-refinement cost proportional to the fixed
query Top-N rather than total library size. This should permit rapid new
co-crystal queries while retaining exact provenance and avoiding guessed
chemistry during SDF/docking export.

## One-time companion contents

- Atomic number, formal charge and aromatic flag for every heavy atom.
- Heavy-heavy bonds with order and aromatic/ring flags.
- Exact heavy-atom membership for every artifact-v1 pharmacophore feature.
- Signed donor/acceptor projection directions and axial aromatic normals when
  geometrically defined, with explicit validity/kind fields.
- Strict non-ring rotatable torsions, their four defining atoms and the bounded
  terminal-side heavy atoms moved by refinement.
- Stable ID and molecule/conformer identity checks against artifact v1.

## Query-time fixed order

1. Load only detailed rigid-refinement Top-N IDs by mmap.
2. Reconstruct chemically identified candidate coordinates and the stored rigid
   pose without reading MOL2.
3. Add projected-color evidence while retaining atom-centered scores.
4. Evaluate element-aware soft vdW penetration against the exact receptor.
5. For at most two smallest terminal torsions, evaluate a locked small angle
   set with bounded beam search; preserve pre/post score, angles and strain
   proxy.
6. Export chemically correct structures only from companion-backed topology.

## Locked initial limits

- Post-refinement candidates: caller-controlled fixed Top-N, never all loose
  pharmacophore hits.
- Maximum terminal torsions optimized per pose: 2.
- Maximum moved heavy atoms per stored terminal torsion: 12.
- Angle offsets: -60, -30, 0, +30 and +60 degrees.
- Beam width: 8.
- Donor/acceptor directions are signed; aromatic normals are axial and compare
  with absolute cosine.
- Protein clash uses element-specific vdW radii and reports raw penetration;
  thresholds remain annotations until enrichment/redocking calibration.

## Confirmatory acceptance criteria

- Companion global IDs and molecule/conformer identities exactly equal artifact
  v1 in original order.
- Atomic numbers, charges, bonds and feature memberships reproduce source RDKit
  molecules on a locked test set.
- Direction vectors are finite unit vectors when valid; invalid directions are
  explicit and contribute no fabricated directional evidence.
- Rigid transforms rotate directions without translation.
- vdW penetration is zero for separated controls and positive with analytically
  expected ordering for controlled overlaps.
- Torsion rotation preserves bond lengths and all nonmoving atom coordinates;
  zero-angle refinement is exactly the input pose and score cannot decrease.
- Repeated companion builds and query refinements are deterministic and hash
  stable.
- Existing v1/FAISS/pharmacophore/Gaussian tests remain unchanged and pass.

## Decision boundary

The new scores are reversible ranking and QC evidence. They are not binding
free energies, experimental affinities or permission to discard the protected
retrieval union. Any hard clash or directional cutoff requires preregistered
redocking/enrichment calibration.
