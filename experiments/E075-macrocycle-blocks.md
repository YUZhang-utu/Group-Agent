# E075: first conformer-level macrocycle blocks

Hypothesis: a conservative ring/amide-state partition followed by capacity-bounded
torsion splits provides reproducible conformer blocks without merging source identities.

Scope: offline pilot, at most 100000 sampled conformers. No production search change,
no score-based block rejection, no building-block inference from source names.

Protocol: read paired coordinate/chemical artifacts and join each identity to the
read-only registry. Detect one unambiguous macrocycle or validate a supplied ordered
heavy-atom ring map. Canonicalize cyclic traversal; compute signed ring torsions,
amide states and eligible single-bond count. Separate hard groups; recursively split
on the largest-variance sin/cos component at the median until capacity is met.
Persist source identities, ring maps, descriptors, membership, split prototypes and
unassigned reasons. Freeze sample seed, parameters and input catalog hashes.

Tests: rigid-transform and atom-reindexing invariance; cis/trans boundaries;
ambiguous graphs; capacity and hard-group separation; deterministic membership;
source identity preservation and invalid-input handling. Record test results before
workstation acceptance. This establishes engineering properties, not enrichment,
energy barriers, balanced chemical diversity or full-library scalability.

## Results (2026-09-25)

Engineering checks: 8 focused tests passed, including existing block exploration.
Real-data extraction: two source prefixes of 1000 records each; 1750 assigned,
250 held unassigned for fused/bridged ring mapping; seven blocks at capacity 1000.
All assigned records have 18 ring atoms, six amides and six peptide units; 584
assigned source molecule names. Output and implementation hashes are recorded in
`data/e075-source-final-2000/report.json`.
Exploratory outcome: bounded extraction and partitioning work on real records;
automatic ring recovery needs broader graph handling. No recall or speedup measured.
