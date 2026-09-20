# E034: expanded-library WEE1 retrieval, Gaussian refinement and E031

Historical implementation handoff. Subsequent workstation results are recorded in
[E034 results](E034_RESULT_2026-09-18.md). At handoff, E033 reported25,813,808
conformers and8,318,351source-grouped molecules;10000candidates/nprobe128 was a
provisional eight-query calibration setting, not E034 acceptance.

Historical local validation:166passed/2dependency skips, including19new E034 tests
and actual mmap Gaussian/E031 integration. Bash syntax passed. The local Windows
interpreter lacked RDKit/FAISS; mock-index checks were not real FAISS validation.

## Run on the workstation

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e034_expanded_wee1.sh
```

Default inputs and output:

- E032: `/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`
- E033: `/mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604`
- QT9: `data/e019_query_8bju/{8BJU.cif,QT9.cif,query_manifest.json}`
- 824: `data/e026_query_1x8b/{1X8B.cif,824.cif,query_manifest.json}`
- Output: `/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/20260918-run1`

These are the existing E026 local query directories; data are not in Git. Override
actual locations without renaming unrelated query files:

```bash
E034_QT9_DIR=/actual/qt9-query-dir \
E034_1X8B_DIR=/actual/1x8b-query-dir \
bash scripts/run_e034_expanded_wee1.sh
```

The runner neither downloads structures nor rebuilds the library. Missing queries,
incomplete library/evidence, failed calibration, identity or lineage mismatches
stop execution. Dependencies match E032/E033/E026: Python3.11+, NumPy, RDKit,
Biopython, Gemmi and FAISS. Defaults:16Gaussian workers,20FAISS threads, single-thread
BLAS; override with `E034_WORKERS`/`E034_THREADS`. Many mmap shards need a larger
file-descriptor limit. The process raises its soft limit within the existing hard
limit, inherited by children; it does not change system configuration.

## Execution and evidence boundaries

1. Validate E033 completion/report/protocol, catalogs, transform/index hashes and
   current registry/artifact/chemical IDs. Historical full payload acceptance comes
   from E033; this rechecks metadata/identity and index bytes, not original MOL2.
2. Copy locked inputs and prepare crystal-specific queries. USRCAT uses standardized
   parent chemistry; Gaussian uses the CCD ligand form from the same crystal instance.
   No verified external-query-to-library molecule-ID mapping means no claim of
   complete chemical self-match removal.
3. Scan exact whole-library USRCAT Top1000 for each query. Test budget10000 at
   nprobe128/256; choose128 if both queries have strict recall>=95%, otherwise256.
   If neither passes, write retrieval evidence and exit2 without refinement.
4. Gaussian-score all selected conformers, refine the union of per-objective
   Top5000 with the existing512pair-seed cap. Do not first collapse to one conformer
   per molecule.
5. E031 annotates three existing rigid objective poses: coverage, per-anchor matches,
   correlation and Top-K overlap. It changes neither ranking nor admission and does
   not perform docking.
6. Report retrieval, fetch, exact-reference, Gaussian and E031 costs/retention
   separately. E031's <=10% gate uses same-query/same-run Gaussian refinement time.
   Execution complete does not imply latency or scientific acceptance.

Target descriptor calibration is not Gaussian or biological recall. Budget reduction
is not chemical rejection. No extra pharmacophore channel is added in this experiment.
The legacy E031 runner no longer gates against the old fixed77.99-second baseline.

## Resume

```bash
bash scripts/run_e034_expanded_wee1.sh --resume
```

A nonexistent output starts fresh even with `--resume`. Existing output without
`protocol.json` is preserved and rejected; use a fresh `E034_OUTPUT` or the real
previous run. Inputs, code, dependencies, host and worker/thread settings must match.
Completed stages retain original timings after hash checks. Interrupted Gaussian
stages reuse validated chunks, but partial refinement reuse makes the latency gate
unavailable. For a new complete timing:

```bash
E034_OUTPUT=/mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/timing-run2 \
bash scripts/run_e034_expanded_wee1.sh
```

Inputs/results are not automatically deleted. A stage without its committed receipt
may repeat after interruption. `FAILED.json` is historical; inspect `RUN_STATUS.json`,
the final marker and report for current status.

## Missing E033 completion marker

`report.md` is written before the completion marker and does not prove successful
completion. A wrong directory, incomplete copy or finalization interruption can
cause a missing marker. Do not fabricate `EVALUATION_COMPLETE.json`. Point
`E033_OUTPUT` at complete evidence, or run:

```bash
bash scripts/run_e034_with_fresh_acceptance.sh
```

This creates timestamp/PID-specific E033/E034 directories, reruns acceptance and
calibration, then starts E034 only on success. It preserves the library and old
reports. Record the printed paths and set `E033_OUTPUT`/`E034_OUTPUT` for resume.

## Return files

Return `report.md` and `report.json`. `retrieval/` holds candidates and exact distances.
Each query has copied inputs, `gaussian/refine/merged-scores.npz`,
`interaction-matches.npz` and manifests. Root `protocol.json` pins inputs/code;
`*.stage.json` are receipts. Startup, integrity and exact-reference scans are not
online query latency; descriptive small-panel timing is not an SLA.
