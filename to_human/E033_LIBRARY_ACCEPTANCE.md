# E033: library acceptance and retrieval calibration

## Workstation execution

Use the Python environment that ran E032:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e033_library_acceptance.sh
```

Default library: `/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`.
Results go to a separate timestamped sibling `e033-library-acceptance-*` directory.
The runner does not rebuild the index, alter source MOL2/registry/precompute data,
or upload library data.

Registration alone does not mean E032 is complete. If `COMPLETE.json` is absent,
check whether feature generation, index training or merging is still running.
Wait for that process, or resume the interrupted original command below. Do not
start two precompute processes at once.

```bash
set -o pipefail
python scripts/precompute_library_batch.py \
  --source /mnt/local/hand/yuzhang/aidd/mc_data \
  --old-artifacts "$PWD/data/e019_artifacts/catalog.json" \
  --old-chemical /mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json \
  --output /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --workers 20 --preserve-index-conflicts --run \
  2>&1 | tee -a /mnt/local/hand/yuzhang/aidd/library-precompute-20260911/run.log
```

## Measurements

1. Library acceptance: completion markers, source registration coverage,
   records = inserted + exact duplicates, counts across stages, contiguous global
   IDs, registry/artifact/chemical identity, manifest lineage, full generated-file
   and pharmacophore-posting SHA256, final index dimensions and count. Historical
   source paths are matched using unique stems and source hashes; no manifest
   rewrite is required. Original MOL2 bytes are not rehashed, as the report states.
2. Timing: index loading, first FAISS call, warm p50/p95, candidate fetch and exact
   USRCAT reranking, full exact-reference scanning, and peak process RSS. Integrity
   checks and reference generation are outside online latency. First-call timing
   is not guaranteed cold-cache timing; small-panel p95 is not a service SLA.
3. Retrieval quality: exact standardized-USRCAT squared L2 across the entire
   library, compared with approximate candidates at Top100 and Top1000. Boundary
   ties are reported separately. A single shard is never treated as global truth.
4. Retention: conformer and source-grouped molecule counts for budgets
   1000/10000/100000 and nprobe64/128/256. Molecule-cap scenarios of100/1000/10000
   retain one descriptor-ranked representative per molecule and report exact
   neighbor molecule coverage. These are capacity scenarios, not scientifically
   accepted docking shortlists or a final single-pose policy.

The default panel deterministically selects eight distinct molecules from uniformly
sampled conformer IDs; molecules with more conformers have greater selection
probability. All conformers of each query molecule are excluded from both exact
truth and retrieval. This is an internal calibration panel, not a biological
holdout. The engineering gate is strict Top1000 recall >=95% for every query.
Select the smallest passing candidate budget and then compare timing. If none
passes, report that explicitly.

## Outputs

- `report.md`: timing, worst-query recall, retained molecules and reduction ranges.
- `report.json`: full parameters, provenance, per-query results and provisional settings.
- `metrics.csv`: query/configuration/cap comparisons.
- `acceptance.json`: integrity acceptance.
- `protocol.json`, `queries.npz`, `truth_*.npz`: frozen settings, queries and exact truth.
- `candidates_*.npz`: global IDs, molecule IDs, exact descriptor distances and cap
  representatives. Uncapped `global_ids` can feed Gaussian only with the same query
  and library; internal-panel candidates cannot be assigned to a WEE1 crystal query.
- `EVALUATION_COMPLETE.json`: execution completed; check `calibration_status`
  separately. `FAILED.json` preserves failure details and partial outputs.

Return `report.md` and `report.json`; an average reduction alone omits worst-query
recall and remaining molecule counts.

## Interpretation

This measures descriptor coarse retrieval, not Gaussian pose quality, interaction
preservation, activity enrichment or docking success. Budget-induced reduction is
not chemical rejection. Next use calibrated settings with actual WEE1 queries and
E024 refinement, measuring its separate cost/retention. Biological enrichment needs
active/decoy labels. Source-name collision stereochemistry QC remains independent.

## Optional settings

```bash
# Use a fresh output and a larger query panel.
E033_THREADS=20 E033_OUTPUT=/mnt/local/hand/yuzhang/aidd/e033-review-2 \
  bash scripts/run_e033_library_acceptance.sh --query-count 16

# Supply independently prepared query descriptors.
bash scripts/run_e033_library_acceptance.sh --queries /absolute/path/queries.npz
```

External NPZ files require raw, non-z-scored `vectors` (n x60) and `names` (n), with
optional `exclude_molecule_ids` (n). Use the same USRCAT generation method and record
provenance; arbitrary 60-dimensional embeddings are not interchangeable.
`--metadata-only` skips large-file hashes and cannot establish full byte acceptance.
The default full check is the formal acceptance path. Both modes scan exact truth.

Method references: [FAISS performance](https://github.com/facebookresearch/faiss/wiki/How-to-make-Faiss-run-faster)
and [FAISS FAQ](https://github.com/facebookresearch/faiss/wiki/FAQ).
The parameters and gate are E033 protocol choices, not vendor guarantees.
