# Benchmarks

## E019 conformer coverage

- Dataset: Platinum high-quality protein-bound ligand conformations.
- Primary metric: symmetry-corrected minimum heavy-atom RMSD from every
  generated ensemble to its bound conformation.
- Report: median, 75th/90th/95th percentiles, CDF at 0.5/1.0/1.5/2.0/2.5
  angstrom, failures, runtime, and conformer count.
- Strata: rotatable-bond count, heavy atoms, molecular weight, ring-system
  count, macrocycle flag, and query-chemotype-matched subset.
- Purpose: calibrate retrieval thresholds and retained top-N; it is not a
  go/no-go test because the production library conformers are immutable.

## E019 retrieval correctness and speed

- Truth: exhaustive joint Gaussian shape/color optimization on a tractable,
  preregistered subset.
- Index metrics: FAISS candidate recall@100/@1000/@10000 versus exhaustive joint
  scoring, candidate reduction, bytes/conformer, and missed-hit analysis across
  nprobe/PQ settings and L1 top-N budgets.
- Latency: build throughput; warm/cold query p50/p95/p99; FAISS, local-record
  fetch, seed/pre-score, local refinement, and torsion-refinement separately.
  Compare random ID-order reads with file-offset-sorted batched reads.
- Scale points: 100K correctness set, 1M prototype, then extrapolation checked
  against at least one larger physical shard before native-engine approval.

## E019 enrichment

- Primary: LIT-PCBA, because experimentally tested inactive compounds reduce
  synthetic-decoy bias.
- Secondary: DUD-E, reported with an explicit bias caveat.
- Metrics: per-target and macro-averaged EF1%, BEDROC (fixed documented alpha),
  hit counts/rates at fixed top 1,000 and 10,000, ROC-AUC only as a secondary
  metric, hit diversity, and confidence intervals. Compare methods relatively;
  LIT-PCBA inactive labels are assay-specific rather than proof of nonbinding.
- Baselines in fixed order: Morgan 2D; USRCAT; aligned shape-only Gaussian;
  shape+unweighted atom-centered color; full anchors+projected color.
- Ablations: anchors, projected polar points, ColorTversky, pocket exclusion,
  and terminal-torsion refinement.

## Decision rules

- Conformer median RMSD <1.0 angstrom: tighter thresholds/smaller top-N may be
  tested. Between 1.0 and 2.0: widen flexible-molecule retention. Above 2.0:
  test lower shape contribution and greater anchor/color reliance, but accept
  only after enrichment validation.
- No native billion-scale implementation until the FAISS/L1.5 route
  demonstrates high exhaustive top-N recall and the additional color/anchor
  machinery improves early enrichment over shape-only.
- Drop anchor weighting/projection if its improvement over unweighted color is
  not significant and robust; this also permits narrower geometric bins.
- Fit thresholds only on development targets and freeze before held-out
  evaluation. Report uncertainty of the random 99.9th percentile; 100,000
  random molecules yield only about 100 expected tail observations.

## Local storage and ingestion budget

- NFS is read-only cold storage and is scanned sequentially once per source
  shard. Online queries never read it.
- Local files are uncompressed: int16 `coords.bin`, int16 `feats.bin`, binary
  `meta.bin`, bounding-box/PMI records, raw float32 USRCAT, and FAISS IVF-PQ.
- Build all per-conformer artifacts during the same shard pass. Atomic done
  markers include source and output checksums and permit restart after reboot.
- Preserve raw USRCAT vectors after FAISS construction so PQ/nlist/nprobe can be
  changed without reparsing NFS.
- Account for FAISS IDs: at one billion vectors, 30-byte PQ codes plus 8-byte
  IDs require about 38 GB before inverted-list/alignment overhead. Budget
  42--46 GB steady-state for m=30 and 52--56 GB with query/L2/Python workspaces,
  leaving roughly 8--12 GB page cache in 64 GB RAM until measured. Standard
  FAISS IVF uses 64-bit `idx_t`; int32 IDs require a custom/implicit layout.
  Assert direct_map is disabled.
- Compare m20 (20-byte codes, contiguous 3-D subquantizers), m30 (30-byte,
  2-D subquantizers), and padded-d64/m32. Benchmark build peak separately from
  a clean-process reload after persistence.

## Result artifact schema

- One immutable directory per query and configuration hash.
- E022 Gaussian batch acceptance additionally requires exact input/output ID
  equality in original order, separate transforms per named objective, raw
  query/candidate/cross overlaps, source hashes, and deterministic numeric
  reruns. Runtime is reported separately for query packaging, artifact reads,
  seed generation, and scoring; synthetic correctness is not an enrichment or
  production-throughput claim.
- E023 performance reporting uses the same fixed candidate IDs for 1- and
  16-worker runs and records wall/user/system time, aggregate CPU utilization,
  peak RSS, chunks reused/computed, candidates/s, and score-array equality.
  A deliberate interruption/restart must reuse completed real chunks before a
  full run is accepted as resumable.
- `query_manifest.json` defines the ordered anchor list using stable IDs, ligand
  atom indices, feature type, atom/projected coordinates, interaction evidence
  (distance, angle, burial), and materialized weights.
- Ranked binary/Parquet rows contain `mol_id`, `conf_id`, and named objective
  poses. Each pose has its own 4x4 transform plus every shape/color score
  evaluated under that transform. Expected objectives include shape-only,
  atom-centered joint, and projected joint; cross-pose score mixing is invalid.
- Each pose stores unnormalized `per_anchor_overlap_raw[]` and
  `per_anchor_self_overlap_raw[]`, mapped by the manifest's ordered stable IDs.
  Anchor fields are reversible reranking views, never admission gates.
- Manifest: source/index hashes, USRCAT mean/std, post-z-score group weights,
  score direction/definitions, calibrated cutoffs, versions, and cache/timing
  conditions.
- Export transformed SDF for top 1,000 only by default.
- L1.5 hard filters are restricted to anchor-independent loose volume/PMI
  bounds. Pharmacophore counts may annotate or reorder but may not discard.
- Default query weights are anchor 1.5, ordinary 1.0, solvent-exposed 0.5.
- Report Spearman rho and top-1,000 overlap between otherwise identical
  unweighted and anchored rankings; inspect discordant hits below 60% overlap.

## E024 scale gates

- Ranked-shard merge correctness: compare heap Top-K byte-for-byte with an
  exhaustive concatenate/sort reference, including tied scores.
- Resource metrics: wall time, peak RSS, page faults and bytes/read at 10M and
  100M synthetic metadata, followed by a physical billion-vector run.
- Retrieval correctness: sweep shard-local Top-K and verify global recall@K
  against frozen exhaustive truth; global merging cannot repair shallow shards.
- Schedule invariants: baseline prefix equality, per-channel new-admission
  budgets, duplicate provenance flags, deterministic output hash.
- Gaussian invariants: E024 slim coarse scores and streaming Top-N must reproduce
  E023 on QT9; refined detailed arrays must reproduce for identical selected IDs.
- Resume: interrupt both coarse and refine stages and require valid chunk reuse
  with no `.partial` files.
