# E110: rank blocks with the same similarity criterion

This benchmark tests whether 500 random molecules per block predict the yield of unseen qualifying molecules under the frozen E109 similarity evaluator. It does not rerun the full library. E094 is the default first experiment; additional schemes are explicit because 500 samples per block can still be expensive.

## Run on the workstation

From the existing Group-Agent checkout, activate the same scientific environment used for E109:

```bash
git pull --ff-only origin feature/structure-guided-chat
export PYTHON=/home/y/yuzhang/envs/dl4s/bin/python
export MDM2_BENCHMARK_OUT="$HOME/mdm2-e109"
export SIMILARITY_BENCHMARK_OUT="$HOME/mdm2-e110-e094"
export HARDWARE="TY2603035.utu.fi; CPU=Intel Core Ultra 9 285; RAM=62Gi"
export WORKERS=20
bash scripts/run_mdm2_similarity_benchmark.sh all
```

The scientific modules are unchanged, preserving E109 prepared code hashes. Do not regenerate E109 preparation. The original workstation prepared and pilot-full directories, including their databases and protocol, are required; the small returned reports archive alone is insufficient.

To inspect cost before starting refinement, run the stages separately:

```bash
bash scripts/run_mdm2_similarity_benchmark.sh plan
# Read new_conformers and estimated_new_refinement_hours in plan/report.json.
bash scripts/run_mdm2_similarity_benchmark.sh fit
bash scripts/run_mdm2_similarity_benchmark.sh replay
```

Use either all or the three stages, not both from the start. Plan/replay require new output directories and refuse overwriting. Fit supports identical resume and reuses verified completed chunks; rerun fit after an interruption, then replay. Failed planning leaves an incomplete directory: use a new output path after correcting the cause. Keep script, query, inputs and worker settings fixed while resuming.

For all three schemes explicitly set SCHEMES="E094 E095 E096" and SIMILARITY_BENCHMARK_OUT="$HOME/mdm2-e110-all" before a new run. This samples 500 per block in every scheme and evaluates the union once; it may require substantially more work. Shared fitting cost is charged conservatively to each similarity arm. Separate scheme runs are preferable when comparing deployment costs.

## Execution and interpretation

1. Plan: verify E109 receipts, sample 500 distinct eligible molecules per block, form the union of their conformers, exclude all sampled molecules from pilot validation, and identify already evaluated conformers (including failures). Planning scans the prepared metadata and may take minutes.
2. Fit: evaluate only missing conformers with the unchanged threshold engine. A sampled molecule is positive in a block if any of its conformers in that block passes stage 7 for any template. Count each molecule once. Rank blocks by positive count, breaking ties by block ID. All conformers of sampled molecules are evaluated for consistent reuse; out-of-block positives do not inflate a block's count.
3. Replay: compare frozen similarity priorities, ChemPLP, EquiScore and 20 random orders on identical remaining pilot molecules. Search cost counts conformers. Candidate yield deduplicates molecules, and fixed budgets include complete blocks only. It does not launch additional refinement.

The main endpoint is qualified candidate yield, not recovery of global RRF Top-K. Default work budgets are 10%, 20%, 50%, 100%; candidate targets are 1000, 10000, 50000. Targets larger than available validation positives are reported as unreached, never assigned an artificial speedup.

Runtime cannot be promised before the actual sample union is counted. The plan estimates new refinement from E109's measured throughput. Fit reports actual elapsed work but includes cache reuse; replay reports conformer work rather than measured search latency. cold_conformer_work adds the complete logical sample cost even when cached. It measures this bounded validation workload, not the cost of a full-library deployment. Sampling positives are excluded from validation yield; this isolates prediction rather than claiming end-to-end production delivery counts. E109 preparation and historical docking-priority costs are separate. A later fresh, paired timing experiment is required before claiming wall-time acceleration.

This is an exploratory follow-up because E109 pilot results motivated it. Sample and validation molecules are disjoint, but the pilot was already inspected at aggregate level. Freeze a successful method before a disjoint confirmatory experiment.

## Return files

The runner writes $SIMILARITY_BENCHMARK_OUT/reports.tar.gz after replay. Return that archive. Main files:

- replay/summary.md: quick comparison at 20% work.
- replay/yield_at_work.csv: candidate yield at fixed budgets, including random repetitions.
- replay/work_to_candidates.csv: work to fixed candidate targets with explicit shortfalls.
- replay/yield_curves.csv: complete prefix curves.
- fit/priorities.csv: learned block priorities and sample success counts.
- plan/report.json and fit/report.json: sample size, reuse, new work and timing provenance.

Protocol: experiments/E110-similarity-aligned-routing.md. Only synthetic tests are run locally; no new scientific result is asserted by this implementation.
