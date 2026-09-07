# E024 Billion-scale tiered retrieval and Gaussian execution protocol

Status: core execution architecture implemented and offline validated

## Motivation

E023 proved that 299,999 conformers can be PCA-scored in about 41 seconds and
that 7,704 selected conformers can be pair-refined in about 78 seconds on the
workstation. A billion-conformer library changes the execution problem: even
cheap per-conformer work and a monolithic result array become expensive, while
the scientifically useful distinction between broad and precise retrieval must
not be lost.

## Hypothesis

A shard-native, fixed-budget hierarchy can preserve broad and precise search
evidence while keeping query-time memory and expensive Gaussian work dependent
on configured budgets rather than total library size. Ranked shard outputs can
be merged exactly with a bounded k-way heap. Coarse Gaussian scores can remain
available as immutable slim chunks, while only the stable union of per-objective
Top-N candidates receives detailed pair-seed refinement.

## Locked invariants

1. FAISS/L1 baseline recall is never reduced by pharmacophore or anchor logic.
2. Pharmacophore strict, balanced, and loose channels only add candidates and
   provenance. Duplicate IDs merge their evidence flags.
3. Each retrieval channel has an explicit admission budget. Expensive work is
   bounded by the sum of these budgets, not by the number of posting-list hits.
4. Full broad evidence is retained by immutable per-shard ranked files, counts,
   hashes, and references. It need not be expanded into one billion-row NPZ.
5. Global Top-K over independently ranked shards is deterministic: higher score
   first, then lower stable global ID for ties.
6. Coarse and refined Gaussian products are both retained. Coarse products use
   slim chunk files containing stable ID and the three float32 objective scores;
   detailed transforms and overlap primitives are recomputed and retained only
   for the refinement selection.
7. Streaming Top-N selection over coarse chunks must equal a full stable sort.
8. Every chunk and run is hash/config locked, atomically promoted, resumable,
   and append-only with respect to older library shards.
9. The query manifest reports both evidence cardinality and compute admission;
   candidate budgeting must never be described as biological rejection.
10. No physical billion-scale latency or recall claim is made without a real
    large-shard benchmark against frozen truth or a declared surrogate.

## Execution hierarchy

1. Each immutable library shard owns its reusable artifact files, FAISS index,
   and pharmacophore postings.
2. Query each FAISS shard for a fixed local Top-K and merge to an exact global
   baseline Top-K without concatenating all shard hits.
3. Produce separately ranked strict, balanced, and loose pharmacophore channel
   files under fixed producer caps. Their complete source manifests remain
   referenced even when online admission is smaller.
4. Build the Gaussian schedule in this order: baseline, strict additions,
   balanced additions, loose exploration additions. Each channel's admitted
   prefix and all duplicate provenance are retained.
5. Run PCA-only Gaussian scoring over the fixed schedule into slim resumable
   chunks. Do not create a monolithic detailed coarse archive.
6. Stream the three objective Top-N sets from the slim chunks and take their
   deterministic union.
7. Run detailed pair-seed Gaussian refinement on that union and preserve the
   existing objective-specific transforms and raw overlap primitives.

## Confirmatory acceptance criteria

- K-way shard merge equals exhaustive concatenation/sort on synthetic shards,
  including score ties, with memory-mapped shard arrays and at most
  `global_K + max_local_K + shard_count` resident working rows.
- Schedule construction preserves the complete admitted FAISS prefix, respects
  every channel budget, merges duplicates, and is deterministic across reruns.
- Anchor/pharmacophore inputs can never remove a baseline ID.
- Slim coarse chunks contain exactly ID plus three named float32 scores and use
  no object arrays; their logical payload is 20 bytes per conformer.
- Streaming per-objective Top-N and selected union equal the existing E023
  monolithic stable-sort result.
- Coarse chunks remain after refinement and detailed refined output contains
  exactly the selected IDs.
- Valid chunks are reused; missing or corrupt chunks alone are recomputed.
- Appending a new ranked shard does not mutate old shard files and changes only
  the new global merge/run manifest.
- Full dependency-light test suite remains green.

## Scale validation ladder

- Offline synthetic correctness and bounded-state tests.
- Existing 299,999-conformer QT9 regression against E023 scores/selections.
- Physical 10M and 100M metadata/index simulations measuring merge latency,
  peak RSS, result bytes/row, and resume behavior.
- At least one true multi-shard production-scale run before enabling a
  `billion_scale_validated` capability flag.
- A physical billion-vector benchmark is required before publishing billion
  throughput, latency, or recall numbers.
