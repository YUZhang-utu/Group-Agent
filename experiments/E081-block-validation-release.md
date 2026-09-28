# E081: independently validate and publish the conformer block workflow

Protocol, 2026-09-28, before new validation runs.

Add a mandatory offline validator to the workstation block wrapper. Exhaustively
check descriptor-to-membership identity and coverage, exported manifest equality,
hard-stratum isolation, capacity, leaf statistics and frozen-tree routing. Validate
the source audit accounting and retain review records as an explicit release hold.
Recompute descriptors and independent plane-normal omega angles for a deterministic
bounded raw-source sample. Export those samples for human inspection. Every fitted
conformer participates in structural checks; sampled source checks are labeled.

Tamper tests must catch incorrect membership even when output hashes are updated,
stale or swapped identities, wrong capacities, mixed strata and changed sources.
Test the real E078 mixed fixture with fresh outputs and all requested capacities.
Do not rerun the E077 full source audit. These tests establish implementation and
data integrity, not biochemical quality or improved retrieval recall.

Publish the full E075-E081 implementation, tests and workstation instructions to
the existing GitHub feature branch, after focused regression, language and diff
checks. Preserve unrelated build output and the earlier E066 patch. Do not upload
raw library sources, large result databases or local credentials.

Results: 56 focused tests passed. Full repository regression: 613 passed, 3 skipped.
Three existing E080 real mixed-fixture builds passed the new validator into fresh
data/e081-final-validation-{backbone,chemistry,typed} directories. Each checked
24 conformers across capacities 4/10000/20000/30000 and independently recomputed
all 24 source geometries. Updated-hash corruption tests rejected wrong leaves,
identity/stratum changes, missing exports, wrong populations, radii and capacity.
These are confirmatory implementation checks and exploratory small-source evidence;
workstation full-library execution and scientific retrieval quality remain pending.
