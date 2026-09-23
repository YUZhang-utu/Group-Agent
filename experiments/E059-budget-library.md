# E059: molecule-budget library retrieval and docking handoff

Protocol, 2026-09-23, before execution. Exploratory, not an activity validation.

Use the existing full-library USRCAT FAISS index with current target crystal
queries. Retrieve unique molecules per template, progressively increasing the
conformer search depth. Fuse retrieval ranks and retain an explicit molecule
budget (default 1,000,000). Load every stored conformer of these molecules and
score against every selected template. No empirical size, feature, extent,
anchor-count, mandatory-contact or composite-score gate in this new mode.
Retain the recorded receptor collision rule. Seed enumeration remains bounded.

Keep best actual pose per molecule/template, annotate dual-region occupancy,
rank within each template, and fuse ranks with RRF (k=60). Keep the complete
ranked pool for later pages; export 100,000 unique molecules per page from the
original, hash-verified MOL2 records, with original names and pose evidence.
The first page is a capacity selection, not a prediction of 100,000 actives.

Checks: zero-contact and below-threshold pose retention; unchanged threshold
mode; molecule deduplication, stable rank ties and pagination; source identity
and hash validation; interrupted chunk resume; explicit insufficient-pool counts.
Compare existing pilot molecule-best scores to crystal self-control scores,
noting selection truncation and different pose search protocols.

ANN recall on MDM2 and whole-run wall time are unknown. WEE1 recall parameters
do not establish either. The workstation run is required for actual data.
