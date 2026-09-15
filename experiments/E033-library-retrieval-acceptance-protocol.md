# E033 — Completed-library acceptance and 3D coarse retrieval calibration

Date: 2026-09-15. Status: protocol before workstation measurement.

Question: after E032, what speed/descriptor-recall/candidate-reduction tradeoff
does the actual expanded library achieve before expensive pose refinement?

Hypothesis: increasing IVF search effort and candidate budget improves coverage
of exhaustive standardized-USRCAT neighbors, with a measurable latency cost.
This is a directed calibration, not an unconstrained parameter sweep.

## Frozen first-run design

- Require E032 COMPLETE.json, private registry, three catalogs, final FAISS.
- Acquire existing batch lock read-only; reject a concurrent writer.
- Full artifact/posting/index SHA256 integrity by default; optional metadata
  mode must be labeled incomplete byte-integrity verification.
- Check source-registration coverage, records=inserted+duplicates, registry and
  artifact counts, contiguous global ranges, artifact/chemical identity arrays,
  manifest lineage, vector shape/finiteness and FAISS dimension/count.
- Seed 20260915; sample uniform global conformer IDs until eight distinct
  molecules are represented (molecule sampling is conformer-count-weighted),
  excluding ALL conformers of the query molecule in both truth and candidates.
  Optional external raw 60D query vectors can replace this calibration panel.
- Reuse frozen mean/std from E032; no training or library rewriting.
- Exact full-corpus squared L2 truth with chunked direct subtraction, stable
  distance/global-ID tie ordering. No concatenated full-library vector matrix.
- nprobe 64/128/256; candidate budgets 1000/10000/100000, clipped for tiny data.
- Truth K 100 and 1000. Report strict ID recall and boundary-tie-aware coverage;
  denominator is actual eligible truth size, never the requested K on tiny data.
- One untimed warmup and three timed single-query searches per configuration.
  Report first-call separately, then warm p50/p95 (small panel descriptive only),
  index-load time, exact-pass wall time, per-query reference compute time,
  candidate fetch/rerank wall time and process peak RSS where available.
- Report actual retained/rejected conformers AND unique source-name-grouped
  molecule IDs; distinguish budget truncation from scientific exclusion.
- Exact descriptor rerank of returned candidates; molecule caps 100/1000/10000
  with one representative conformer each. Report unique-molecule truth coverage
  at every cap. These are capacity scenarios, not final docking nominations.
- Proposed engineering gate: every panel query >=0.95 strict recall of exact
  top-1000 in candidates. Report smallest tested budget passing, then fastest
  nprobe at that budget; no passing setting means no recommendation.
- Preserve each configuration's raw candidate IDs/distances and capped
  molecule representatives, with query/index/catalog hashes and parameters.

## Scope and interpretation

USRCAT neighbor agreement is not pose-alignment quality, activity enrichment,
stereochemical QC or docking success. Sampled library queries are calibration,
not held-out scientific validation. Gaussian/interaction-filter retention and
WEE1 labeled enrichment remain separate follow-up measurements using E024/E031.
No "typical filtering percentage" is assumed: report measured distributions
for this corpus, panel and explicit budgets. No original library data are sent
to GitHub; only code/protocol/tests are published.

References: FAISS official speed/accuracy parameters and benchmarking caveats:
https://github.com/facebookresearch/faiss/wiki/How-to-make-Faiss-run-faster
https://github.com/facebookresearch/faiss/wiki/FAQ
