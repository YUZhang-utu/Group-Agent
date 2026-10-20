# E115 quick native-state review

Three implemented diagnostics: detailed polar exposure and whole-color
attribution for MOL-1EC0C7A6D657; reference-only spatial contact grouping with
one-to-one assignments; feasibility-preserving local rigid refinement.
Production thresholds are unchanged. This is an outcome-selected 65-pose review,
not another 10000-conformer search and not an estimate of library hit rate.

Local validation: 23 focused tests passed, including grouping/assignment, kernel
replay, feasible-pose retention, polar exposure and export fixtures. CLI help and
English-content guard passed. Real workstation execution and scientific results
are pending; synthetic tests do not establish improvement in retrieval.

Use the already activated aidd_workstation environment and repository root.
The receipt must be the original E113 directory containing pose_diagnostics.json,
report.json and protocol.json, not the E114 output or the compressed archive.

```bash
git pull --ff-only origin feature/structure-guided-chat
python scripts/review_contact_regions.py \
  --previous "$HOME/mdm2-e111-strict-10000" \
  --prepared "$HOME/mdm2-e109/prepared" \
  --reference "$HOME/mdm2-e109/pilot-full" \
  --receipt "$E113_RECEIPT" \
  --output "$HOME/mdm2-e115-contact-review" \
  --seconds 180
tar -czf "$HOME/e115-contact-review-report.tar.gz" \
  -C "$HOME/mdm2-e115-contact-review" .
```

Set E113_RECEIPT to the same directory used successfully for E114. Use a new output
directory. Optimization has a cooperative three-minute budget; source hashes,
native preparation, mandatory endpoint replay and primary-candidate exposure add
overhead. Actual runtime requires a workstation measurement.

Read report.json for completion and optimization budget flags. reviewed_poses.json
preserves both endpoints, every feasible trial, and a selected feasible pose.
Its grouped score includes per-group assignments; inspect missing polar and
aromatic groups rather than interpreting the average alone. identity_controls.json
checks supported native self-matches, while leave-self-out rows provide the
cross-reference comparison. A grouped mean is not an affinity or acceptance label.

primary_review.json contains per-feature typed overlap, feature-to-atom membership,
isolated/bound polar SASA, and distance-only complementary polar-neighbor hypotheses.
Indices are zero-based; SDF atom serials are index+1. These neighbors are not
angle-validated hydrogen bonds. Polar burial alone does not establish satisfaction,
permeability or binding. Receptor scope is the verified native target chain only.
SDF exports preserve companion bonds and charges, but do not reconstruct stereo tags.
Run view_primary.pml from the output directory in PyMOL to inspect the selected pose.

The next decision requires returned data: determine whether grouped scoring fixes
redundant points without hiding missing essential chemistry, and whether improved
shape retains contact and geometry constraints. Do not relax color or launch the
full library solely because this diagnostic produces higher scores.
