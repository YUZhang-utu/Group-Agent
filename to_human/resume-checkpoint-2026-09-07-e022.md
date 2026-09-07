# Resume checkpoint — E022 artifact-backed Gaussian reranking

## Verified state

- Production library: 100,000 molecules / 299,999 conformers.
- Immutable conformer artifacts and incremental FAISS use stable global IDs.
- QT9 incremental FAISS recall: top-1,000 100%; top-10,000 99.69%.
- E020 partial pharmacophore retrieval is implemented with nested
  loose/balanced/strict tiers and lossless FAISS union.
- E021 Gaussian overlap/Tversky/rigid-seed mathematical kernel is validated.
- E022 now reads candidate chemistry directly from artifact mmaps, builds a
  hashed co-crystal query package, and stores independent shape-only,
  atom-centered unweighted joint, and anchor-weighted joint best poses.
- Anchors only rerank and never remove an L1 candidate.
- Full dependency-light test suite: 90 passed at implementation commit
  `05cb2f5`.

## Next workstation actions

1. Finish/verify transfer of reusable artifacts and databases to local NVMe.
2. Regenerate the artifact catalog on Linux so its shard paths are local.
3. Run the real E020 QT9 validator and retain `pharmacophore-hits.npz`.
4. Run `prepare-gaussian-query` for 8BJU/QT9.
5. Run the locked first 1,000-candidate Gaussian runtime test, inspect outputs,
   then decide the full L1 execution size from measured throughput.
6. Record real runtime/ranking/enrichment before adding projected color,
   exclusion volume, or terminal torsion.

Commands and output definitions are in
`docs/gaussian-artifact-reranking.md`.

## Scientific boundary

Artifact schema v1 has atom-centered features only. Projected candidate color,
atomic radii/topology, exclusion volume, and torsional refinement are pending;
no enrichment or production-throughput claim has been made. Keep unrelated
necessity/enhancement work outside this project and its evidence trail.
