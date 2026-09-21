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

## E025 multi-query aggregation gates

- Conformer collapse: compare every objective-specific winner with an exhaustive
  group-by molecule reference, including tied scores.
- Query union: measure unique-query hits, multi-query support distribution and
  overlap of protected-query versus consensus admissions.
- Rank robustness: apply monotonic score rescaling independently per query and
  require unchanged within-site consensus order.
- Site isolation: identical molecule IDs in different pockets remain separate
  admission keys and never contribute to cross-site RRF.
- Docking expansion: every task references one admitted molecule, one immutable
  query result, one exact receptor and one retained conformer/transform.

## E031 key-interaction sidecar gates

- Correctness: compare deterministic assignment with known non-greedy optima;
  reject type, direction-kind, unit-vector, identity, and input-hash mismatch.
- Immutability: hash the rigid result before and after scoring and require exact
  equality; retain sidecar IDs in identical order.
- Speed: measure both WEE1 queries with `/usr/bin/time -v`; compare directional
  wall time with the accepted 77.99-second QT9 rigid-refine baseline only on
  identical hardware/IDs. Target <=10%, with absolute time and page faults.
- Exploratory kernel measurement: 23,112 three-anchor/twelve-feature
  assignments took 3.72 seconds locally (6,217/s), excluding artifact and
  companion mmap reads. This is not a production acceptance measurement.
- Effect: report score min/median/max, median matched-anchor count, Spearman rho
  versus each Gaussian objective, and Top-100/500/1000 overlap. Do not activate
  reranking without held-out multi-target enrichment or redocking evidence.

## E033 expanded-library acceptance
Protocol: experiments/E033-library-retrieval-acceptance-protocol.md. Full generated payload checks, registry/artifact/chemical IDs, source accounting and frozen index lineage precede retrieval evaluation. Eight deterministic distinct-molecule queries exclude all same-molecule conformers. Full-corpus streaming exact USRCAT L2 truth; nprobe 64/128/256 and candidate budgets 1k/10k/100k; strict and boundary-tie coverage of top100/top1000; first-call and warm search latency, fetch/rerank cost, peak RSS, retained conformers and unique molecule IDs. Explicit molecule-cap scenarios retain one descriptor-ranked representative and report reference-molecule coverage. Proposed .95 minimum per-query recall gate is engineering calibration only, not biological or Gaussian pose quality. Report actual budget-induced reduction distributions; no assumed typical rejection percentage.

## E033 reported workstation acceptance — 2026-09-18
This supersedes earlier pending-result descriptions. User-supplied report: acceptance passed on 25,813,808 conformers / 8,318,351 source-grouped molecules. Eight-query calibration at budget 10000, nprobe 128: worst/mean recall .99/.997125; search p50/p95 6.356/7.209 ms; search plus fetch and exact descriptor reranking median 45.945 ms; retained molecules 8352–9834. Provisional engineering setting only: budget-induced reduction does not measure chemical rejection, activity enrichment, or pose quality; panel is not independent holdout. Full raw report/hashes have not been independently reviewed locally. Workstation output: /mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604. Continue using target-matched WEE1 queries and expanded-library catalogs, preserving conformer provenance. Gaussian and E031 costs/retention remain separately unmeasured on this candidate set. See to_human/resume-checkpoint-2026-09-18-e033-panel-passed.md.

## E034 regression and workstation gates — 2026-09-18
Local:166 passed,2 dependency skips; Bash syntax passed. New19 checks cover E033 hash/lineage and insufficient acceptance, stage mutation/resume, both-query .95 gate/fallback/failure, real Gaussian/E031 fixture integration and immutable IDs, partial-resume timing invalidation, output containment, and many-shard file-descriptor limits. Workstation: fresh locked crystal queries; budget10000 and nprobe128/256; exact standardized-USRCAT Top100/1000; separate index/exact/retrieval/coarse/refine/E031 timings; refined counts and retention; E031 annotation only. Cold-cache, biological enrichment and docking success remain outside acceptance. See experiments/E034-expanded-wee1-refinement-protocol.md.

## 2026-09-18 — E034 workstation result received
User supplied completed report.md and partial report.json from /mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118. Both locked WEE1 queries have exact-USRCAT Top1000 recall1.0 at budget10000/nprobe128 and256; select128. QT9 retains6299conformers/2212molecules after refine, refine31.242s, E0314.870s, overhead15.59% fails10% gate.824 retains6791/2961, refine69.954s, E0314.421s, overhead6.32% passes. This is protocol-aligned user-reported engineering evidence, not an independent artifact/hash audit. Exploratory ranking disagreement: anchored Gaussian versus E031 rho .1801/.1919; conformer-level Top100 intersections1/2. Do not promote E031 to ranking or infer biological/pose superiority. Next reuse these artifacts for molecule-level disagreement/pose review and numerical-equivalence-preserving E031 profiling; do not rerun registration or change the gate. Details: to_human/E034_RESULT_2026-09-18.md and .json. No new computation in this result-recording turn.

## E035 local evidence and boundary — 2026-09-18
Reference Hungarian assignment is unchanged. Batched geometry with feature-only reads is opt-in;100000synthetic kernel cases have zero assignment differences, max numerical error1.11e-16 and exactfloat32 rankings. Kernel speedup2.10–3.68x on nonempty groups is promising but cannot substitute for measured workstation end-to-end latency. E035 now exports molecule-level discordance and both winning conformers without changing admission/ranking.188local checks pass/2dependency skips. Real-library100k/1M distinct-conformer constructed-pose equivalence plus the original39270refined-pose score comparisons remain workstation tasks. No biological or pose-quality conclusion follows from kernel equivalence. See to_human/E035_WORKSTATION_RUN.md.

## E036 evaluation — 2026-09-18
Use identical input hashes/samples for serial vs parallel stress comparisons; no reused chunks for wall-speedup measurements. Verify max error<=1e-6, exact assignments and both conformer/molecule ranks. E036 Gaussian scheduling requires all saved arrays exactly equal to E034 (including poses/transforms); E031 remains sidecar-only. Fresh index candidates must exactly reproduce locked schedules. Report startup/index load separately, actual query execution separately, and final equivalence checks separately. Historical E034 timing comparisons are exploratory; workstation speedup pending. Local194passed/2dependency skips, including spawn/partial resume/tamper and full fake-index pipeline. Handoff:to_human/E036_FAST_3D_RUN.md.

## E037 benchmarks — 2026-09-18
Microbenchmark:scripts/benchmark_bounded_seeds.py, seeded20260918,12cases, feature counts8/16/32,5queryanchors,3alternating repeats,cap512. Exact object prefix required; to_human/E037_LOCAL_SEED_BENCHMARK.json records all timings/environment/source hash. Local205passed/2dependency skips. Workstation suite compares4fixed combinations, two reverse-order repeats, same-machine query medians vs old-reference; partial reuse/profiles excluded. Original E034 scores/poses/identities remain exact acceptance reference. No post-hoc budget tuning, no biological claims. Actual end-to-end results pending.

## E038 interface benchmark — 2026-09-19
222passed/2dependency skips, final16targeted checks passed; standalone offline smoke labeled synthetic. Coverage: strict planner schema, auth error redaction, unknown action/path/dependency rejection, Project isolation, UniProt identity/species/ambiguity, AF3 construct/CCD/provenance, shell-free fake scientific runners, artifact tamper and resume. Real LLM plan fidelity, real protein endpoints, AF3 inference and actual GPU timings remain unmeasured. to_human/E038_LOCAL_VALIDATION.json records evidence boundaries. Do not equate offline smoke with live workflow or scientific acceptance.


## E039 software acceptance (2026-09-20)

Run `python scripts/check_english.py` for maintained UTF-8 content and `python -m
pytest -q` with `PYTHONPATH=src` for regression. Observed: no Han text, 224 passed,
two dependency skips. Four repository skills passed the skill-creator validator.
Tests enforce exact skill/action coverage, valid reference links, routing and
English-plan rejection. This is not a live LLM, AF3 or search latency benchmark.


## E040 software acceptance (2026-09-21)

Full regression: 230 passed, two dependency skips. Container tests verify the
supplied mount/argument semantics, prerequisites, failure preservation, reuse and
image provenance. Provider tests verify explicit selection and credential isolation.
This establishes software behavior under fixtures, not API reliability, AF3 quality
or GPU runtime. English-content and shell syntax checks passed.


## E041 conversational acceptance (2026-09-21)

Full regression: 240 passed, two optional dependency skips. Ten chat tests passed
again after final bounded-log/result-display changes. Tests cover persistent
conversation clarification, owned task isolation, report-derived answers and timing,
real child-process blocking, cancellation, restart recovery, existing-run attachment,
HTTP authentication/origin handling and malformed router output. Existing prompt
adapter test observes atomic running/completed progress while AF3 is simulated.
JavaScript and shell syntax, English content and updated skill metadata passed.
Desktop UI rendered in headless Edge and was visually inspected; no live provider
or GPU inference was invoked by these tests. Local preview uses isolated fixtures.
The initial sandboxed Edge renderer failed; approved headless execution produced
the inspected screenshot under outputs/e041/chat-desktop.png.

E042: 258 local tests passed, 2 optional-dependency skips; 15 evidence/selection tests passed after final fixture portability adjustment. Fixture evidence only; live model routing, library exports and manual review pending.

## E043 engineering benchmarks

Exhaustive retrieval is compared to dense exact Top-K on a multi-shard fixture,
including ties and final-shard candidates. Classified-feature lineage and a
mocked reference-pose gate are tested separately. Passing local tests does not
establish real 25,813,808-conformer timing, category accuracy, or licensed docking
compatibility. Exhaustive candidate sets are not required to equal old ANN sets.

## E044 full-library conditions

Coverage requires all contiguous catalog IDs plus accepted distinct molecule counts.
Counts are deduplicated across chunks; all conditions must hold on one conformer
pose before aggregation. Partial runs fail full-coverage acceptance. Local tests
compare the worker's Gaussian transforms against the existing scoring path and
exercise corruption/restart, last-chunk hits and export provenance. No full-library
latency, storage or biological-quality number has been measured for this mode.

## E045 necessary-condition validation

Analytic basis: for positive spatial-times-angular feature score >= t with
angular <= 1, positional error <= sigma*sqrt(-2 log(t)), bounded by cutoff.
Pair-distance discrepancy is bounded by the sum of errors. Tests retain passing
rigid fixtures under random transforms and boundary radii, reject impossible
assignment graphs, and compare final pass IDs with an unfiltered small fixture.
Legacy passing-pose contradictions stop a run. Workstation prefilter retention,
wall-time speedup and scalability have not yet been measured.
