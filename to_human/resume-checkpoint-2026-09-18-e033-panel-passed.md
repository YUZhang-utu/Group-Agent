# E033 calibration passed: handoff to target-specific refinement

Historical checkpoint,2026-09-18. Based on user-pasted report.md and partial JSON,
not local re-execution or independent inspection of full artifacts/hashes/tables.
This superseded the September15waiting-for-registration state.

## Reported results

- Acceptance passed:25813808conformers,8318351source-grouped molecules.
- Eight-query calibration panel, not an independent holdout.
- Provisional budget10000/nprobe128: worst recall0.99, mean0.997125; passes0.95gate.
- Search p50/p95:0.00635649/0.00720884s; search+fetch+exact descriptor rerank median0.04594543s.
- Retained molecules8352-9834; budget-induced reduction99.8818%-99.8996%, not chemical rejection.
- Budget10000/nprobe256: worst0.998, mean0.9995; search p50/p95 0.01076154/0.01169474s,
  search+fetch/rerank median0.05118022s. Retain as a tradeoff comparison, not proof
  of target-specific superiority.
- All1000budgets failed;10000/100000with nprobe64also failed.
- Linux5.14.0-611.49.2.el9_7x86_64,24CPUs,FAISS1.14.3. Evaluation wall231.3171s,
  peakRSS1208819712bytes; neither is per-query online cost.

## Boundaries and handoff

This measures standardized USRCAT recall, not enrichment, Gaussian poses or docking.
First-call timing may be warm; eight-query p95 is descriptive. Source grouping is
not demonstrated full-library chemical deduplication. Collision stereochemistry
QC remains separate.

Confirmed report.md directory:
`/mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604`.
Other raw artifacts had not been inspected locally.

1. Inspect report/acceptance/protocol/metrics and lineage; preserve the accepted library.
2. Retrieve fresh candidates for the actual8BJU/QT9and1X8B/824crystal queries using
   the frozen USRCAT transform; start10000/128and compare10000/256.
3. E033 panel candidates belong only to their original queries. Use actual WEE1
   `global_ids` with matching Gaussian queries. Molecule caps are capacity analysis.
4. Verify expanded artifact/chemical catalogs, IDs, query/receptor identities and
   hashes. Historical E031 defaults pointed at the old library.
5. Measure Gaussian coarse/refine and E031 separately, with retention counts;
   keep E031 annotation-only.
6. Replace the old77.99sbaseline from7704candidates with same-query/same-run timing.

Subsequent implementation note: E034 removed the fixed77.99s gate and added the
new runner; see [E034 handoff](E034_WORKSTATION_RUN.md). Workstation execution was
still pending at this checkpoint. Later results are recorded separately in E034
and E035. Machine-readable E033 summary: `e033-user-reported-20260918.json`.
Resume from the passed25.8M-conformer panel, not the older registration-wait state.
