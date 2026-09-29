# E089: source-backed candidate/reference review

Protocol before execution. User workstation E088 angular upper-bound coverage
is 97.76% at 35 degrees, 99.51% at 45. No property builds found at sibling paths.
Test whether candidate backbone/property/steric discrepancies are comparable to
variation within their definite target classes. Do not infer compatibility from
angle alone or set a cutoff to meet coverage.

Select up to 100 original boundary classes per band (<=35, 35-45, >45 control),
one minimum-hash CID per class/band. Include all proposed rotations and up to
three distinct-molecule reference conformers per target class. References are
indexed first distinct molecules, not unbiased representatives. No minimum-score
cherry-picking. Report every pair and within-class reference controls.

Scan only needed MOL2 files once, parse and recompute chemistry only for selected
records. Verify source hash/name, descriptor reconstruction and residue mapping.
Align new property profiles to stored unit atom indices explicitly, because typed
and backbone cyclic canonicalization may differ. Rotate boundary into reference
order; use proper-rotation Kabsch backbone RMSD and separate property/steric
metrics. Export selected MOL2, provenance, profiles, all pairs and summary.

No automatic reassignment, docking or affinity claims. Whole-library acceptance
and rare-state coverage remain unproven. Failed runs have no complete receipt.

## Local execution result

Seven focused tests passed. Final small source-backed run selected 6 conformers
(3 boundary, 3 reference) and reproduced source descriptors. It emitted 5 boundary
pairs and 1 within-class pair; 1 target lacks a within-class control. All outputs
in data/e089-final-pair-review; this is exploratory fixture evidence only.
