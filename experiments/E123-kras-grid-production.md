# E123 retained-pool grid prioritization and docking export

Protocol recorded before implementation and execution. Evaluate only E121 pool
members, using unchanged E122 receptors (7RPZ/12FF), contacts, three conformers,
eight exact poses and profiles. Do not interpret failed finite searches as
molecular rejection. Use SQLite committed target/molecule results for resume,
with source/code/settings seals. Share conformer reads for overlapping targets.
No coordinate search outside the pool; the compact catalog ID arrays are read
once to map up to three existing conformers per pool molecule.

Export at most 1,000,000 balanced-profile molecules per target ordered by the
supporting balanced pose score, with deterministic molecule-ID ties. This is a
computational quota, not an affinity ranking or a chemical-diversity selection.
Use the corresponding original MOL2 conformer, preserving its coordinates and
bonds. Check source hash, parsing, atom order, coordinate agreement, same-pose
contacts and exact grid profile on original coordinates before priority export.
Failed source checks are quarantined, not activity negatives. Do not silently
backfill the quota. Retain every result and an ID/tier manifest; original E121
exports remain available for all reserve/unresolved IDs. Split at 100,000 records.
Evaluation resumes at committed batches. Export resumes at completed targets;
an interrupted target is restarted in a new attempt directory, never appended.

Tests: pool-only mapping; checkpoint skip and seal mismatch; deterministic quota;
native fixture end-to-end run, source recheck and MOL2 counts; original E122 tests.
Workstation evidence pending. User-reported pilot: 10,000 per target, 97.35 s
compute, approximately 9.08 h full-pool extrapolation excluding export. The two
native controls do not establish recall, enrichment, or candidate pose validity.
