# E025 resume checkpoint - 2026-09-07

## Accepted state

- E024 shard-native fixed-budget search and slim/detailed Gaussian execution are
  implemented; physical 10M/100M/1B scale validation remains pending.
- E025 multi-cocrystal aggregation is implemented with site-scoped union,
  objective-specific conformer-to-molecule collapse, Top-M alternatives, RRF,
  protected query quotas, molecule admission, and receptor-specific tasks.
- Dependency-light suite: 102 passed.
- User ran the complete E025 validator on the real QT9 detailed refinement
  result and reported successful completion with no failed checks.
- The real validation covers the single-query path only. It must not be reported
  as a completed multi-cocrystal RRF validation.

## Next action

Run a second independent co-crystal through E024, add it to the same biological
`site_id` in the plan, execute `aggregate-multi-cocrystal-search`, and record:

- unique molecules per query;
- overlap and query-support distribution;
- protected-only, consensus-only, and overlapping admissions;
- number of unique admitted molecules and receptor-specific docking tasks;
- output/source hash validation using the accepted complete validator.

After that acceptance, begin the versioned docking adapter with co-crystal
redocking as its first protocol. Do not mix the separate KRAS necessity or
enhancement project into this AIDD research state.

