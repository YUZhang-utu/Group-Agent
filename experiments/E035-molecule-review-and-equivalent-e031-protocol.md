# E035 — Molecule-level disagreement, review poses and equivalent E031 acceleration

Freeze before implementation checks, 2026-09-18. User requests all three tasks
plus large-scale validation. Existing E034 outputs remain immutable and E031
remains annotation-only; do not change admissions or Gaussian ranking.

Hypothesis: moving repeated feature transforms/validation/geometry into bounded
NumPy batches and reading only interaction-required fields reduces E031 time
without changing its one-to-one assignment objective or chosen assignments.
Keep the legacy implementation and Hungarian solver as the reference.

Input: completed E034 recheck-20260918-163847-316118, both query packages, rigid
NPZs and E031 sidecars; bind report/completion, query, result and catalog hashes.
Existing 6,299 + 6,791 conformers, three objectives = 39,270 real refined poses.

Molecule review: independently maximize anchored Gaussian and anchored E031
over each source-grouped molecule's conformers. Keep both winning conformer IDs,
cross-scores and exact objective transforms. Stable tie order by global ID.
Compute molecule-level Spearman and Top100/500/1000 overlap. Export up to10
molecules in each group: Gaussian top100 with E031 outside top500; E031 top100
with Gaussian outside top500; both top100. Do not fill missing groups with
weaker examples. Export both representatives when distinct, native-frame
heavy-atom SDF with bonds/charges, anchor assignments/contributions, copied
receptor mmCIF and a PyMOL review script. Human inspection remains pending.

Equivalence gates: IDs/order/objectives/anchor assignments exact; all score and
per-anchor arrays finite, max absolute deviation <=1e-6 with rtol=0; stable
conformer and molecule rankings exact, including ties. Any mismatch fails the
gate and prevents claiming equivalence. Never silently fall back to fewer poses.

Timing: profile one reference call separately (not used as latency); perform
three alternating reference/optimized full calls on each original rigid NPZ.
Include feature I/O, validation, scoring, analysis and output serialization.
Report distributions and median, speedup and ratio to the recorded same-query
E034 refine baseline (QT9 31.24205s;824 69.95412s). This uses a historical
baseline on a verified matching machine, not a new Gaussian run. If machine
identity differs, ratio/gate unavailable. Preserve <=10% gate; no result
claimed until actual workstation measurements. Repeated timings are descriptive.

Scale validation: deterministic stratified sampling across the entire artifact
catalog, default100,000 distinct library conformers, plus optional1,000,000.
Evaluate both WEE1 anchor sets and three deterministic centroid-aligned rigid
orientations per conformer:600,000 comparisons at default scale;6,000,000 at1M.
These are real library features but constructed stress poses, NOT Gaussian
refined poses or docking/enrichment validation. Reference and optimized engines
must agree on every score/assignment; record sampled IDs, per-shard coverage,
feature/anchor counts, mismatch counts, max deviations and throughput. Chunked
resume must preserve hashes; do not repeat old candidates to pretend to scale.
Actual expanded-library execution occurs on the user's workstation; local
randomized and mmap-fixture checks are reported separately.

Offline validation: empty/zero-feature candidates, typed signed/axial direction
semantics, exact/cutoff/tied/near-tied assignments, invalid transforms/identities,
scalar/batch equivalence, differing per-molecule representatives, source
immutability, export coordinate/bond consistency, deterministic unique sampling,
and interruption/resume input/output tamper rejection. Record speed results
even if the optimization fails the original latency goal.
