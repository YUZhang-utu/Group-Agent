# E108: run the label-free routing benchmark

## What this can establish

The scientific claim is **reference-candidate recovery at lower search cost**,
not biological activity, binding affinity, or superiority on every target. The
library has no experimental labels. Three evidence levels are kept separate:

| Level | Available inputs | Output and permissible interpretation |
|---|---|---|
| Historical comparison | Existing million-conformer reference and twelve arms | Preservation of the historical ANN-derived candidates; no new search |
| Controlled routing | Complete evaluated regular universe plus frozen priorities | Reference recovery versus examined conformers, against random ordering |
| Actual deployment | Matched fresh executions and independently checked recovery | Measured time ratio at matched recovery, output budget and conditions |

The first level is ready to run on existing workstation files. The other two
require additional evidence; the software does not invent it or silently launch
searches. The historical reference is **not** exhaustive same-policy scoring.
Do not state that the library screens better merely because the scripts complete.

## 1. Run the existing-data comparison now

Python 3.9+ standard library only; no RDKit, docking, models, NumPy, or installation.
Run from your existing checkout on the workstation:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
bash scripts/run_mdm2_routing_benchmark.sh "$HOME/mdm2-benchmark-e108"
```

Use a new output directory. Never remove existing scientific outputs to rerun the
analysis. The command reads the twelve completed E106 runs from the October 9
handoff and the historical export at:

```text
/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers
```

The original input search was `PROMPT-59addb9ffd6444c3` with eleven templates.
The historical export contains 1,000,000 conformers / 460,332 distinct molecules
according to the supplied receipt. Counts are validated from files, not hardcoded.
The 3,854,146 input conformers were associated with the historical selected million
molecules. This is an ANN-derived candidate subset, not exhaustive pose scoring.

Override locations if files have moved:

```bash
MDM2_REFERENCE=/path/to/mdm2-11templates-top1m-conformers \
MDM2_RUNS_ROOT=/path/to/runs \
bash scripts/run_mdm2_routing_benchmark.sh "$HOME/mdm2-benchmark-e108-new" \
  --registry /path/to/the/original/registry.sqlite3
```

The registry path is normally discovered from the first arm's sealed protocol.
It is opened read-only. Keep the historical `report.json`, `selected-conformers.csv`
and current per-arm `report.json`, `protocol.json`, `query/report.json`,
`selected/report.json`, `selected/candidates.jsonl` accessible. MOL2 files are not
read for this analysis. By default also retain each arm's `members.sqlite` and
`ranking.sqlite` for stage-loss accounting. If these were not transferred, append
`--skip-stages` explicitly; this omits that analysis, not the main comparison.
Do not move just the MOL2 files and lose their mappings.
Expect several GB of available RAM for million-record identity tables; no runtime
promise is made. Avoid a concurrent writer to the source registry or result files.

### Identity safeguard

Current structural `cid` is not necessarily the historical artifact identifier.
The tool joins current `artifact_cid` to historical `conformer_id`, verifies the
original record SHA-256 against the registry, resolves canonical molecule IDs,
and cross-checks shared global IDs. Names, formulae and exported `conf_*` titles
are not identity keys. This avoids false zero overlap and namespace mismatches.

### Files to inspect

| File | Meaning |
|---|---|
| `report.md`, `audit.json` | Scope, provenance, warnings and missing evidence |
| `arm_summary.csv` | Work, survivors, exports, molecules, shortfalls, cache information |
| `historical_recovery.csv` | Reference Top-K recovery at each output prefix, for conformers and molecules |
| `pairwise_overlap.csv` | Exact overlap among all twelve arms, in canonical namespaces |
| `top5_top10_gains.csv` | Added/lost candidates and extra work when block limit increases |
| `reference_stage_coverage.csv` | Reference candidates in selected blocks, visited blocks, qualified pool and final export |
| `recovery_work.svg` | Open directly in a browser; historical recovery versus searched work |

For `historical_recovery.csv`:

- `reference_recall = reference_hits / reference_k`.
- `selected_reference_fraction` is reference membership, **not active precision**.
- Molecule reference ranks are first occurrence in the historical conformer
  ranking, not the older independently constructed molecule ranking.
- `full_output=True` rows describe actual delivered outputs. Prefix rows describe
  sensitivity at a smaller export budget, but retain the full arm's search cost.
- K larger than available reference units is omitted, never silently clipped.
- The six Top10-limit outputs are the primary full-budget comparison. All twelve
  outputs are also compared at their common export prefix (normally 76,933).
- The two Top5 shortfalls are preserved. Do not fill them from other blocks.
- Four Top5/Top10 pairs have identical reported counts; exact equality is tested,
  not inferred. Top10 is a maximum limit, not necessarily ten searched blocks.
- Reference Top-1,000,000 recall has a 10% ceiling for a 100,000-conformer export;
  focus also on Top-1,000/10,000/100,000. Do not interpret this ceiling as failure.

Historical tiered quota/RRF and current threshold-funnel RRF differ. Current RRF
is recomputed across the searched arm's poses, not a universal library rank.
Outside-reference candidates are unknown, not bad candidates. A low overlap can
reflect policy differences, not just routing loss. Raw RRF values are not directly
comparable across candidate universes. No reported cached time is used as speedup.

Stage accounting separates references outside selected blocks, references in blocks
not visited after early stopping, references searched but not qualified under the
current policy, and qualified references dropped by output ranking. This is more
informative than attributing every lost reference to partitioning. Stage accounting
uses exact conformers; molecule-level final recovery is reported separately.

## 2. Stronger controlled benchmark: routing discovery versus work

Predeclare K, a 90% reference-recovery target, independent query panels and budgets
before examining new outcomes. For a feasible first experiment, freeze a
molecule-sampled sublibrary with all its selected conformers. Keep the original
block memberships. Sample across the eligible library, not the union of selected
arms. Reserve molecules used in block-priority fitting; do not score them as held-
out successes. State the exclusion rule and hashes. If historical priorities are
reused without exclusion, label the experiment retrospective/in-sample instead.

Evaluate **every member of this common sublibrary**, including rejected members,
using the same query, thresholds, geometry and reference-ranking policy. Save
the results once; all routing replays reuse them. A sublibrary result does not
prove full-library recall. Full regular-library replay is possible with sufficient
memory; the loader retains one record per universe member and three membership
indexes. Start with the bounded sublibrary rather than loading 25 million rows
on an unprovisioned workstation. Random orders are generated one at a time.

The input schema is deliberately explicit, not a guess at a completed full-library
artifact. No such matched artifact has yet been verified in the local records.

`universe.csv` (every evaluated conformer, original IDs and memberships):

```csv
conformer_id,molecule_id,qualified,reference_rank,E094,E095,E096
artifact-a,molecule-a,1,1,parent-a,property-a,spatial-a
artifact-b,molecule-b,0,,parent-b,property-b,spatial-b
artifact-c,molecule-c,1,2,parent-b,property-b,spatial-b
```

`qualified` is 0 or 1 under the common computational filter, **not activity**.
All qualified conformers must have unique contiguous ranks; rejected conformers
have blank rank. Freeze one rank policy, e.g. the current conformer positive-
contact/RRF policy on the complete evaluated universe. Do not use missing scores
as rejections. Apply deterministic ties. Exclude the special-exhaustive pool from
this regular-block experiment and separately report its size and treatment.

`priorities.csv` contains all six full priority lists, not just leading blocks:

```csv
scheme,scorer,rank,block_id
E094,chemplp,1,parent-a
E094,chemplp,2,parent-b
E094,equiscore,1,parent-b
E094,equiscore,2,parent-a
```

Also supply E095/E096 lists. Ranks must be contiguous per scheme/scorer. Blocks
absent from a sublibrary are allowed and skipped. Every present block must have
a priority. Freeze priorities from training molecules before evaluation; preserve
the existing best-ten-distinct-molecule statistic. Do not refit from test outcomes.

`universe-manifest.json`:

```json
{
  "scope": "complete_fixed_sublibrary",
  "conformers": 3,
  "library_identity": "REPLACE_WITH_SOURCE_AND_SUBSAMPLE_HASHES",
  "query_identity": "REPLACE_WITH_ADOPTED_QUERY_HASH",
  "evaluation_policy": "REPLACE_WITH_FROZEN_ENGINE_THRESHOLDS_RANKING_CODE_HASHES",
  "priority_fit_exclusion": "REPLACE_WITH_MOLECULE_DISJOINT_SPLIT_RULE_AND_HASH",
  "universe_sha256": "REPLACE_WITH_SHA256",
  "priorities_sha256": "REPLACE_WITH_SHA256"
}
```

Get file hashes with `sha256sum universe.csv priorities.csv`. These manifest
statements are provenance declarations, not automatically proven population
completeness. Keep evaluator receipts and the source membership manifest for audit.

```bash
python scripts/benchmark_routing.py replay \
  --universe universe.csv --priorities priorities.csv --manifest universe-manifest.json \
  --top-k 100 1000 10000 --target-recall 0.9 --repeats 20 --seed 20261014 \
  --output "$HOME/mdm2-routing-replay-e108"
```

Outputs:

- `routing_curves.csv`: six trained orders, random block orders for each scheme,
  and a separately labeled uniform-conformer random baseline.
- `recall_at_work.csv`: recovery under 1%, 5%, 10%, 20%, 50%, 100% work budgets.
  Only completed blocks within budget count; no fractional block interpolation.
- `work_to_recall.csv`: actual whole-block work needed to reach the fixed recovery.
- `replay_summary.csv`: median/min/max work across random repeats; these ranges
  describe order variability, not confidence intervals across independent targets.

Compare schemes with the same scorer, scorers within the same scheme, and each
against its random block order. Block count is not a fair compute budget when sizes
differ. Random uniform candidates and random blocks are intentionally different.
Discovery recall here is **before final RRF/export truncation**. It isolates routing
quality, while the historical analysis evaluates delivered outputs. Do not claim
this replay reproduces the current Top5/Top10 deployed export or measures runtime.

## 3. Measure real speed, not cache reuse

The existing Top10 timings reuse results. They cannot be used as fresh speedups.
For an isolated routing comparison, the unblocked baseline and routed method must
share library/subsample, query, eligibility, scoring/rank policy, output budget,
hardware, worker counts and I/O timing scope. Historical ANN versus current
threshold search compares two pipelines, not only the effect of partitioning.

Run complete no-block and routed executions with fresh scientific output locations,
no shared result-cache hits or resume, and the same OS-cache policy. Alternate
execution order over at least three repeats. A warm filesystem cache may be retained
for every method if explicitly stated. Do not claim disk-cold performance. The tool
does not clear caches or stop other jobs. Record competing workload externally.

The explicit `measure` subcommand is the **only command in this toolkit that launches
another program**. Historical/replay/summary commands are analysis only. Provide
the already-approved scientific CLI command as an argv list; this document does not
invent a full-library command or submit new workstation searches automatically.

```json
{
  "method": "E096-chemplp",
  "comparison_id": "mdm2-heldout-panel-v1",
  "hardware": "REPLACE_WITH_CPU_GPU_RAM_HOST",
  "workers": 8,
  "protocol_identity": "REPLACE_WITH_COMMON_LIBRARY_QUERY_ENGINE_BUDGET_HASHES",
  "phase": "online",
  "cache_policy": "fresh_application_cache",
  "repeat": 1,
  "cwd": "/absolute/path/to/Group-Agent",
  "command": ["python", "PATH_TO_APPROVED_RUNNER.py", "--output", "/fresh/scientific/output"],
  "result_receipt": "/fresh/scientific/output/report.json"
}
```

The scientific receipt must explicitly record `status=complete`, `cache_reused=false`
and `resumed_execution=false`; absence is not interpreted as a cold successful run.
For a baseline whose receipt lacks these fields, add them in an audited runner
based on actual execution, not by manually relabeling an old receipt. A prior
existing result receipt is rejected. Command failures retain logs and ineligible
timing evidence instead of silently disappearing.

```bash
python scripts/benchmark_routing.py measure --spec run-spec.json --output /fresh/timing-output
python scripts/benchmark_routing.py timing-summary --baseline no-block \
  --receipts /path/to/timing-runs/*/audit.json --output /fresh/timing-summary
```

`timing-summary` requires at least three identically labeled repeats per method and
identical declared conditions; it reports median paired baseline/method time ratios.
It does **not** verify that a faster method achieved the required recall/output
budget. Join the recovery result from the same scientific output before interpreting
the ratio. At matched work, compare recovery; at matched recovery and output cap,
compare runtime. Never compare a 76,933-output shortfall as a completed 100,000 run.

Measure one-time partition/index building and priority sampling/scoring as `setup`,
separately from `online`. Do not pool phases in timing summaries. Report:

```text
total cost for Q queries = setup cost + sum(online query costs)
amortized cost per query = total cost / Q
```

More blocks may require more initial scoring. Preserve that cost. Do not compare
kernel time or summed worker CPU time with end-to-end wall time.

## Decision and next evidence to return

First return the historical `audit.json`, CSVs and SVG, not all MOL2 files.
They identify candidate-preservation/work tradeoffs and any identity problems.
Then freeze the matched universe experiment and collect timing receipts if a
stronger screening-efficiency claim is needed. Extend to additional independent
queries/targets before claiming general superiority. No activity labels are needed
for this retrieval-efficiency claim, but the complete reference and fair costs are.

Avoid declaring a single winner from high overlap alone. Choose a predeclared
recall/output budget, inspect the Pareto tradeoff between retained candidates,
work, measured runtime and molecule coverage, and report failures/shortfalls.
Initial protocol: `experiments/E108-label-free-routing-benchmark.md`.
