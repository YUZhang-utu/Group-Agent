# E034 — Expanded-library WEE1 retrieval, Gaussian refinement and E031

Protocol date: 2026-09-18. Freeze this protocol in Git before implementation validation.

Question: does the E033 provisional 10,000-candidate / nprobe=128 setting retain
at least 95% of exact standardized-USRCAT Top-1000 for each of the locked
8BJU:QT9:A:601 and 1X8B:824:A:901 queries, and what are the separate costs and
retained counts of Gaussian refinement and annotation-only E031?

Inputs: the E032 completed batch and the E033 report directory
`/mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604`.
Verify E033 completion/report checksum, full acceptance, passing calibration,
catalog and transform hashes, current index checksum and catalog/chemical
lineage. Keep the accepted library read-only under its existing batch lock.
Reuse local locked query manifests and mmCIF/CCD files; do not fetch or silently
substitute a query. Snapshot query inputs and record hashes in a new output.

Retrieval: prepare actual crystal-coordinate USRCAT using the existing parent
standardization routine; build Gaussian packages from the same source and exact
ligand instance. Scan the expanded descriptor corpus once for both exact
Top-1000 references. Compare nprobe=128 and 256 at budget 10000, with repeats=5,
Top-100 and Top-1000 strict/tie-aware recall, timings and unique molecule counts.
External WEE1 queries do not receive library-ID self-exclusion unless a verified
mapping exists; record this explicitly. Do not reuse E033 library-panel candidates.

Selection: choose 128 if both queries pass .95 Top-1000 strict recall; otherwise
choose 256 only if both pass; if neither passes, save results and stop before
Gaussian refinement. Do not silently raise budgets or relax the gate.

Refinement: keep all retrieved conformers (no molecule representative cap);
existing scaled Gaussian coarse stage then union of Top-5000 per objective,
max_pair_seeds=512, coarse chunks=2000, refine chunks=250, default 16 workers.
Run E031 on each resulting immutable rigid pose file, preserve all three
objective-specific poses and record analysis without changing ranking.

Measure: counts and molecule retention at retrieval/refinement; FAISS query,
fetch/exact-rerank, exact reference scan, coarse, refine and E031 costs separately.
Compare full E031 call wall time against this run's same-query Gaussian refine
stage wall time, not the old 77.99 s baseline. <=10% is an engineering target;
a failure is a valid reported result, not reason to alter scores. Do not evaluate
the timing gate on partially reused Gaussian work; report unavailable instead.

Resume: explicit reuse of the same output, frozen parameters/input hashes and
per-stage output hashes; reuse completed stages, delegate interrupted Gaussian
chunks to its existing validated resume. Reject altered inputs or completed
outputs. Keep original timings with reused completed stages. No implicit deletion.

Offline checks: lineage/report tampering, changed queries/output, candidate
identity/range, recall failure gating, correct baseline/timing-resume policy,
real Gaussian/E031 fixture integration, deterministic stage reuse. Real 25.8M
run stays pending for the user's Linux workstation. Report tests separately
from real retrieval/pose/biological validation.

No activity enrichment, docking acceptance, validated chemistry rejection rate
or production latency guarantee is claimed by this experiment.
