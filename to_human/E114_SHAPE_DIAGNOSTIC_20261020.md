# E114: why the shape floor rejects the saved candidates

## Conclusion

The E113 zero-hit observation is reproducible for the examined original seeds, but does not establish absence of shape-compatible molecules among the 10000 conformers. A bounded local rigid refinement reveals an alignment-search limitation: three library candidate poses and four cross-reference control poses newly exceed the unchanged shape 0.50 floor. These are NOT accepted hits; receptor contacts, color and clashes after refinement remain pending workstation recheck.

## Evidence and protocol

Protocol commit b674b05 precedes analysis. E113 receipt hashes verified. Nine reference coordinate clouds reconstructed from recovered self-controls with shape near one. Recomputed original shapes agree within 1e-6 for all 65 panel poses, supporting coordinate/formula consistency. Panel comprises 35 geometry-qualified near misses, top 20 library shape poses and ten leave-self-out best-shape poses. One leave-self-out ligand had no contact-qualified saved pose. This is an outcome-selected diagnostic panel; it cannot estimate population recall.

Final local run: setup 1.836 s, optimization 1.267 s, total 3.103 s. One bounded L-BFGS-B rigid refinement per pose, <=180 objective calls, rotation-vector components +/-0.5 radians, translation components +/-2 A. No torsion changes, docking or global guarantee. Retain original pose if no improvement. Some optimizations end at the evaluation cap; endpoints are best observed, not asserted global optima.

## Results

| Panel | Poses | Newly cross shape 0.50 | Highest refined shape |
|---|---:|---:|---:|
| Original geometry-qualified near misses | 35 | 0 | 0.484843 |
| Top original library shape poses | 20 | 3 | 0.513249 |
| Leave-self-out crystal controls | 10 | 4 | 0.974120 |

Library crossings:

| Molecule | Reference | Before | After |
|---|---|---:|---:|
| MOL-31E0CA0A163F | 7BJ6 | 0.458384 | 0.504565 |
| MOL-1EC0C7A6D657 | 6Q9L | 0.438796 | 0.513249 |
| MOL-627753FFDDA8 | 6Q9L | 0.428254 | 0.507699 |

Control crossings: 6Q9L against 7BJ6 0.329411 -> 0.601383; 5C5A against 7BMG 0.372866 -> 0.565268; 5LN2 against 7BJ6 0.424966 -> 0.529662; 7NA2 against 4OBA 0.488599 -> 0.574207. They show local shape-search headroom even for known ligands, not guaranteed compatible binding in the alternate receptor.

## Size and representation

All 35 original geometry-qualified candidates are larger in heavy-atom count than their matched reference: ratios 1.13 to 1.70. For the six previously highlighted coverage>=0.5 candidates, ratios are approximately 1.27 to 1.57. Their original poses place 63-73% of reference atom centers within 1.5 A of a candidate atom, but only 41-49% of candidate atom centers within 1.5 A of a reference atom. This establishes asymmetric whole-ligand overlap. It does not locate those extra atoms relative to solvent or establish that they are harmless.

Example MOL-D9DD9F80898F: 57 versus 41 heavy atoms, original shape 0.358446, locally refined 0.430628. MOL-1637D1A5CD6B: 58 versus 37, 0.288312 -> 0.321202. None of the six crosses 0.50. Thus their original rejection is not fully rescued by this small rigid optimization.

For every panel case, the self-overlap-norm Cauchy-Schwarz upper bound on Tanimoto is above 0.50. There is no demonstrated mathematical size-only ceiling below 0.50; poor placement/geometry also matters. A loose upper bound does not prove an attainable alignment.

Kernel-width sensitivity is substantial: MOL-D9DD9F80898F at the unchanged original pose scores 0.240 at sigma 0.75 A, 0.358 at 1.0 A, and 0.527 at 1.5 A. This is smoothing sensitivity, not evidence to adopt sigma 1.5. Thresholds cannot be transferred across kernel definitions or to external shape-volume implementations without calibration.

## Decision

Keep sigma 1.0 and shape 0.50 for the next check. Do not conclude that the standard is wrong or that 10000 molecules contain no potential binders. The known limitation is finite seed alignment without local refinement; the unresolved question is whether improved shape is compatible with required contacts and receptor geometry.

Only 65 original/refined pose pairs need workstation endpoint evaluation (55 library proposals, ten controls), with no new library search. The supplied checker verifies E113 lineage, original pose transforms and replayed original/refined shape scores, then recalculates contacts, color, collision and coverage on the same endpoint. Its counts are pose diagnostics, not full-library rates.

## Workstation check

Use the activated path-based aidd_workstation environment and repository root. Pull feature/structure-guided-chat. Set E113_RECEIPT to the existing E113 output directory containing report.json and pose_diagnostics.json; do not pass the tar file itself.

    python scripts/check_refined_shape_poses.py --previous "$HOME/mdm2-e111-strict-10000" --prepared "$HOME/mdm2-e109/prepared" --reference "$HOME/mdm2-e109/pilot-full" --receipt "$E113_RECEIPT" --output "$HOME/mdm2-e114-shape-check" --seconds 120

The default proposals are to_human/E114_SHAPE_PROPOSALS.json. Use a new output directory. Evaluation has a cooperative 120-second cap; loading, hashes and native setup add overhead. Return:

    tar -czf "$HOME/e114-shape-check-report.tar.gz" -C "$HOME/mdm2-e114-shape-check" .

## Artifacts

Local diagnostics: data/e114-shape-diagnostic-final-20261020. Saved proposal bundle: to_human/E114_SHAPE_PROPOSALS.json. Scripts: scripts/diagnose_saved_shapes.py and scripts/check_refined_shape_poses.py. Local focused tests cover Gaussian norm bounds, rigid invariance, transform composition, refinement recovery, asymmetric coverage, timeout handling and endpoint rejection of a shape-perfect but clashing pose.
