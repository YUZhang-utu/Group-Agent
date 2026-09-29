# E091: calibrated candidate-specific boundary routing

Protocol before implementation/run. Seek a small exhaustive special set, with
1% as a provisional engineering target, not a classification threshold. Never
tune compatibility to force coverage. Preserve original identities and known
cis/trans/chirality compatibility from E088. Route <=35 and 35-45 separately;
>45 remains special in this revision. No whole-library recall claim.

Collect up to 64 distinct-molecule reference conformers per candidate target,
using a deterministic hash of molecule ID (not first database entries). Split
each class into disjoint fit/calibration/check sets 50/25/25%. Require at least
32 references. Select up to 12 diverse fit prototypes by farthest-first on the
continuous property, steric and backbone torsion features, keeping calibration
and check molecules out. Sparse classes remain ineligible for confident routing.

Calibrate a joint nearest-prototype score on calibration molecules: scaled
property RMS, steric RMS, and circular backbone sin/cos RMS, with a maximum
across channels and a single supporting prototype. Use a fixed 95th percentile
of calibration scores. Report independent check retention; require >=80% for
operational routing. This check is used as a gate, not an unbiased final test.
Future independent search evaluation is still needed.

Stream selected raw MOL2 records, verify identity and original descriptors,
and compute compact double-precision features for all <=45 candidates and
references. No need to recompute 24M conformers, but ~1.26M boundary records
plus references require substantial chemistry computation. Record exact source
hashes, candidate rotation, supporting prototype, normalized score and all
compatible alternatives. Unsupported cases stay special with explicit reasons.
Frozen source artifacts remain unchanged; new files form a sidecar assignment.

Verify conservation of every boundary CID, stable definite membership, exact
per-band counts, source evidence, disjoint molecular splits, jointly supported
acceptance, proper rotation handling, and unchanged old files. Whole-library
result/coverage remains pending workstation execution.

## Confirmed target and local evidence

User confirmed <=1% as the engineering target during this turn. At the admitted
library size this permits 246,637 special conformers. Eleven focused tests passed;
real small-data execution preserves source SQLite hashes and all identities.
Its sparse references correctly leave all three boundary records special; no
forced coverage. New full-library assignment remains pending workstation run.
