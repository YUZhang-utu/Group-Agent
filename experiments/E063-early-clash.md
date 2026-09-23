# E063: exact early clash rejection and explicit worker budgets

2026-09-23. User requests early pocket checks and 20 workers with 4048-conformer
chunks. A centroid or pocket-sphere membership rule is not equivalent to atom
clash exclusion and is not adopted. For zero allowed clashes only, test eight
actual heavy atoms first, in bounded batches; reject a pose as soon as any probe
clashes. Survivors undergo the complete original predicate. Positive allowed
clash fractions use the complete batched check. No seeds or candidates removed
without a proven physical failure. Seed generation still runs unchanged.

Compare masks against the scalar reference across randomized transforms,
threshold boundaries, hard exclusions and nonzero clash fractions. Verify fewer
KD query points on a deliberately clashing panel, and exact Chat worker/chunk
parameter propagation. Run full regression. Real speed and CPU scaling require
workstation measurement; no linear 16-to-20-worker speedup claim.

Result: 475 passed, 2 skipped; scalar-mask equality, exact strict boundary,
hard exclusions, nonzero fractions and explicit workers/chunk propagation pass.
Constructed 10x60-atom test: 80 probe points plus 60 survivor points, versus
600 full-query points. Work count only; workstation speedup is unmeasured.
