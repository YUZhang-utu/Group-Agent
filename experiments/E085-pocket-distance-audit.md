# E085: pocket distance audit (2026-09-28)

Protocol recorded before recalculation. Exploratory analysis, not a calibrated
target-independent classifier. Preserve E084 outputs unchanged.

Hypotheses: four inactive fields result from atom-centred radius 2 A being
inside the VDW-plus-probe exclusion; local zeros are censored by the 30 A^3
trigger; the old C1 may contain meaningful substructure but its apparent modes
must be checked after removing censoring. PDB support is not occupancy.

Verify 22IZ against original search cache, CIF identity and official RCSB.
Exclude the entire 7BJ0 entry by explicit user instruction, preserving files.
Use the same alignment, seed and admitted structures for a controlled comparison.
Define nondirectional fields as free-space contact shells extending feature_radius
beyond each atom's VDW-plus-probe exclusion. Preserve donor/acceptor direction
proxies. Report source feature counts and all channel volumes, not just distances.

Retain local raw maximum IoU loss, legacy hard threshold and smooth support
version: max(region IoU loss * changed_volume/(changed_volume+30 A^3)).
The smooth version is exploratory; audit all three instead of assuming zero
means identical. Sweep cluster cutoff 0.20 to 0.60 in 0.05 steps; examine shell
width 1, 2, 3 A and local support scale 10, 30, 50 A^3 independently.
Record complete memberships, representative radius violations and C6 membership.
Inspect old C1 member strata and local-distance subclusters; do not infer two
physical states merely from a thresholded one-dimensional distribution.
No automatic C6 deletion or adoption, no parameter selection to recover six.

Validation: analytic feature-support regression for all six channels; local
subthreshold change regression; whole-entry exclusion; existing pocket/workflow
tests; checksum raw results and deliver independent CSV tables and matrices.

## E085 completed exploratory result

67 chains/59 entries retained; 7BJ0 excluded entry-wide. Geometry arrays match
E084 exactly for retained chains. Confirmed four-field support bug; all 67 now
have hydrophobic/aromatic support. C1 zero/nonzero strata 43/12 are weakly
separated; uncensored zeros reduce to the representative itself. Corrected
chemical + legacy local cutoffs .30/.40/.50 give 9/4/2; smooth gives 3/1/1.
72 settings exported. No final state count, no C6 removal, no user adoption.
See to_human/E085_POCKET_DISTANCE_AUDIT.md and E085_VALIDATION.json.
