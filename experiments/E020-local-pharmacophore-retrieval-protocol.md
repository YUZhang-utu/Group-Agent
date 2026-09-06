# E020 local pharmacophore retrieval protocol

Status: implementation completed -- offline confirmatory validation passed;
real QT9 workstation benchmark pending

## Hypothesis

A shard-local invariant index of pharmacophore feature-type pairs and binned
distances can be computed once from the existing conformer artifacts and reused
for every co-crystal query. Taking its candidates as a union with the existing
anchor-independent FAISS pool should recover partial three-dimensional motifs
without rescanning MOL2 or allowing anchors to reduce baseline recall.

## Protocol

1. Read only immutable `meta.bin` and `feats.bin` shard artifacts.
2. For every conformer, enumerate unique unordered compatible feature-type
   pairs up to a configured maximum distance and encode their distance bins.
3. Persist sorted `(pair_key, global_id)` postings, a key/offset table, hashes,
   source-manifest hash, feature definitions, bin width and distance limit.
4. Index every shard independently. Adding a library shard must not rewrite
   previous pharmacophore shard indices.
5. Compile each query manifest once into typed anchor pairs and distance-bin
   windows for loose, balanced and strict profiles.
6. Query the immutable postings without reading MOL2. Record which query pairs
   admitted each conformer and assign nested evidence tiers.
7. Union every local-pharmacophore result with externally supplied L1/FAISS
   IDs. No profile is permitted to remove an external L1 candidate.
8. Persist a query-specific result manifest with input/index hashes and exact
   profile definitions.

## Prespecified profiles

- `loose`: +/- 2.0 angstrom pair tolerance; one matched query pair.
- `balanced`: +/- 1.0 angstrom; up to two matched query pairs.
- `strict`: +/- 0.5 angstrom; up to three matched query pairs.

When a query contains fewer pairs, the required match count is capped by the
number available. A one-anchor query has no invariant 3D geometry and must fall
back to other broad channels rather than fabricate spatial evidence.

## Acceptance criteria

- Repeated builds are byte-stable for the same input.
- Query tiers are nested: strict is a subset of balanced, balanced of loose.
- A synthetic rotated/translated candidate is recovered from invariant pair
  geometry.
- Partial one-pair matches appear in loose even when they fail stricter tiers.
- External L1 IDs survive unchanged in the final union.
- Source/index/query hashes and matching-pair evidence are serialized.
- Existing tests remain green.

## Classification

The implementation tests are confirmatory against this protocol. Performance
and recall on the complete QT9/library artifacts remain a separate workstation
experiment and must not be inferred from unit tests.

## Offline results (2026-09-06)

- Implemented independent immutable pair-postings indices per artifact shard.
- Implemented query compilation without library access and hashed query plans.
- Implemented loose/balanced/strict nested retrieval and an anchor-safe union
  with externally supplied FAISS/L1 IDs.
- A synthetic translated query recovered the exact invariant three-pair match;
  a one-pair partial match entered loose only; an unrelated external L1 ID was
  retained unchanged.
- Complete dependency-light suite: 78 passed.
- The full 299,999-conformer index build and QT9 recall/runtime measurements
  remain unexecuted on the Linux workstation and are not claimed here.
