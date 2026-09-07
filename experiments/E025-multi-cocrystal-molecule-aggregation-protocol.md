# E025 Multi-cocrystal 3D-search and molecule aggregation protocol

Status: implemented and offline validated

## Offline result (2026-09-07)

- Implemented site-scoped query union, independent objective-specific
  conformer collapse, Top-M alternative conformers, within-query molecule
  ranks, rank percentiles, RRF consensus, protected query lanes, and unique
  molecule docking admission.
- Implemented receptor/query-specific docking-task export carrying the selected
  conformer and candidate-to-query transform.
- Synthetic tests confirmed unique-query hits survive, multi-query support is
  counted, objective-specific conformers remain independent, monotonic raw-score
  rescaling cannot change fusion order, sites remain separate, checksums are
  enforced, and identical reruns reproduce output hashes.
- Result is an offline orchestration validation. Real multi-cocrystal campaign
  aggregation and docking-engine execution remain workstation tasks.

## Motivation

The pre-docking search path currently produces provenance-rich results for one
co-crystal query at a time. A target can have several useful co-crystals that
represent different ligand chemotypes and receptor conformations. Intersecting
their hits would discard query-specific opportunities, while averaging raw
Gaussian scores would mix non-comparable score distributions.

## Hypothesis

Independent query runs can be combined losslessly at conformer, molecule, and
docking-task levels by using stable molecule identity, within-query ranks, and
site-scoped reciprocal-rank fusion. Per-query protected quotas plus a consensus
lane should retain both conformation-specific and repeatedly supported hits.

## Locked input contract

A versioned JSON plan lists every independent query with:

- unique `query_id`;
- `receptor_id` identifying the exact docking receptor;
- `site_id` grouping only biologically equivalent pockets;
- path to one completed detailed Gaussian result NPZ;
- optional result-manifest path and descriptive metadata.

The detailed result must contain stable global, molecule and conformer IDs,
objective names, per-objective scores and 4x4 transforms. Every input file and
optional manifest is hashed. Duplicate query IDs fail.

## Locked aggregation rules

1. Query runs remain immutable and independently auditable.
2. Different `site_id` groups are never fused into one score or admission list.
3. Within each query and objective, conformers collapse to one best conformer
   per molecule by score descending, then stable global ID ascending. The top-M
   primary-objective conformers remain attached as alternatives.
4. The configurable primary objective defaults to
   `atomcentered_anchored_joint`; all other objective-specific best poses remain
   in evidence and are not overwritten.
5. Across queries in the same site, aggregation is a union. Missing support
   never deletes a molecule.
6. Raw scores are not averaged across queries. Molecules are compared using
   one-based within-query molecule ranks, rank percentiles, and reciprocal rank
   fusion `sum(1 / (rrf_k + rank))`.
7. `query_support_count` and `receptor_support_count` are evidence fields, not
   hard filters.
8. Docking admission first takes a protected Top-N lane from every query, then
   fills from site-level RRF consensus up to a per-site global limit.
9. The per-site global limit must be at least
   `query_count * per_query_quota`; otherwise query protection is impossible
   and the run fails rather than silently violating it.
10. Every admitted molecule creates one docking task per supporting query,
    carrying receptor, site, selected conformer, source result, objective,
    score, rank, and candidate-to-query transform.

## Outputs

- `query-molecule-evidence.jsonl`: one row per query/molecule with all
  objective-specific best poses and Top-M primary conformers;
- `molecule-summary.jsonl`: one row per site/molecule with support, best rank,
  best percentile, RRF and source queries;
- `docking-admission.jsonl`: unique admitted molecules and reversible reasons;
- `docking-tasks.jsonl`: receptor/query-specific docking jobs;
- `manifest.json`: configuration, hashes, counts and invariant checks.

## Confirmatory acceptance criteria

- A molecule found by only one query survives the site union.
- Multiple conformers collapse to the correct best conformer with deterministic
  tie handling, while Top-M alternatives remain present.
- Different objectives retain independent best conformers/transforms.
- Raw score rescaling in one query does not change rank-fusion ordering when
  its within-query molecule order is unchanged.
- Different sites produce independent summaries and admission budgets.
- Every query receives its protected quota when enough molecules are present.
- Admission contains unique site/molecule keys; docking tasks preserve every
  supporting receptor/query and reference an admitted molecule.
- Repeated identical runs have identical data-file hashes.
- Invalid, duplicate, non-finite, mismatched or non-detailed inputs fail.
- Full dependency-light test suite remains green.
