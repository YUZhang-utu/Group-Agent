# E044: Uncapped full-library condition counts

User requirement: inspect every stored conformer against the existing classified
query features, retain match evidence and count source-grouped molecules. Neither
descriptor Top-K nor Gaussian Top-N may restrict membership. Existing diagnostic
thresholds 0.25, 0.5 and 0.75 remain exploratory display levels, not acceptance
cutoffs. No additional threshold choice is required before computing the table.

Hypothesis: bounded batches with persistent workers can evaluate all catalog IDs
without global candidate arrays or global score merges. Chunk receipts and
disk-backed molecule aggregation allow interruption without double counting.

Protocol: use a sealed classification report as the query/evidence definition;
validate contiguous catalog ranges and shard integrity; generate the existing
bounded rigid-pose seeds for every conformer, retaining the anchored Gaussian
objective pose; score all classified features on that pose. Persist scores,
assignments, transforms and identities by chunk. Count each feature/threshold
using source-grouped molecule IDs, merging duplicate conformers across chunks.
All/any conditions must hold on the same saved conformer pose before molecule
deduplication. Feature aliases share a score column and are not independent tests.

No descriptor ranking rejection, Gaussian Top-N, molecule cap or silent failures.
Partial/pilot output is never a full-library count. Failure stops completion.
This enumerates all library conformers, not all possible orientations or torsions;
pose generation remains heuristic and can miss a feature-compatible pose.

Validation before release: synthetic last-chunk hits, duplicate molecules across
chunks, missing/gapped IDs, same-pose conditions, resumed counting, corrupted
receipts, partial coverage and worker scoring parity. Record real-library timing,
storage and chemistry validation as pending until workstation execution.
