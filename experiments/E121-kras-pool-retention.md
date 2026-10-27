# E121: retained direct-contact pool handoff

## Protocol (before execution)

Hypothesis: E120 local rigid-search failures do not justify molecule-level
rejection. Reuse all completed direct-pool IDs and committed E120 successes,
partitioning each target into priority and unresolved candidates without loss.

Inputs: completed direct counts.sqlite and its protocol/report; an optional
E120 ranking.sqlite snapshot from a stopped job; the sealed source registry.
No full-library pose search or geometry-array traversal is required for this
handoff. Resolve one original conformer per pool ID by indexed registry lookup.
Use the saved best conformer for E120-supported IDs; otherwise use the lowest
registry conformer_index (a deterministic representative, not a scored pose).

All unresolved IDs remain retained. There is deliberately no unvalidated new
chemical rejection rule. Source parse failures are quarantined, not declared
nonbinders. Failed exact-source pose rechecks downgrade priority to unresolved.
No molecule quota; MOL2 batches contain at most 100,000 records.

Validation: fixture union equals input pool for each target; tiers are disjoint;
outside-pool ranking IDs fail; source checksums and identity mismatch fail;
strict source validation and downgrade preserve accounting; snapshot includes
committed transactions even when original report still says running. Preserve
raw source records except title. Tests establish software behavior, not recall.

Time scope: mapping and raw-file export I/O only. Several million MOL2 records
can be large and source files may require sequential traversal. Do not promise
same-day completion without actual throughput. New flexibility-aware screening
is a separate experiment, not silently substituted into this handoff.
