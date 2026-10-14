# E109: matched block-routing benchmark

## What this answers

Compare three partitions (E094/E095/E096) and two existing priority scorers
(ChemPLP/EquiScore). Top5/10/20 are operating points, not separately trained
methods. A single complete evaluation supplies the common reference for all six
methods and repeated random block orders. Only visited IDs change; query,
qualification, pose generation and ranking stay fixed.

This is label-free recovery of computational candidates, not experimental activity
prediction. Do not label reference Top-K as known actives, claim biological
precision, or assume speedup improves with library size before measurements.

The adopted query includes eleven co-crystal templates, Gaussian shape/feature
color, optional protein-contact matching and the 5C5A pocket clash gate. Existing
7NA2_A block priorities remain frozen. These are different receptor contexts and
are explicitly recorded, not silently presented as an ensemble. No new water,
metal, spatial-group or occupancy constraints are enabled. Old ANN Top-million
exports are not the matched reference. The current threshold evaluator is reused.

## Workstation setup

Use the existing AIDD environment with RDKit, NumPy and its screening dependencies.
No local synthetic test substitutes for a real evaluator pilot. Source library,
sampling databases, imported rankings and the completed current campaign must
remain available at their original paths. The default campaign is the recorded
PROMPT-5f7c28e4bd704616 block campaign; override `MDM2_CAMPAIGN` if relocated.
Invalid/missing seals fail rather than silently rebuilding a different benchmark.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
export MDM2_BENCHMARK_OUT="$HOME/mdm2-e109"
export HARDWARE="replace with actual CPU, RAM, storage and host description"
export WORKERS=20
bash scripts/run_mdm2_matched_benchmark.sh prepare
bash scripts/run_mdm2_matched_benchmark.sh pilot
```

Preparation maps the full common regular population, verifies membership counts
and saves six complete priority lists. All conformers of molecules used for block
priority fitting are excluded by default. Special/review populations are out of
scope. Mapping uses the existing artifact mapper, which holds large metadata sets
in RAM; allow ample memory and disk. This stage is not a cheap Top20 search.

The default pilot evaluates a deterministic 1% molecule subset twice: full and
all E094 blocks. It requires identical evaluated IDs, best-pose evidence and ordered
candidate CSVs, then replays all six strategies. Inspect `pilot-equivalence/report.json`
and `pilot-replay/`. Smaller pilots can use `FRACTION=0.001`, but must retain at
least 100 qualified candidates for the default replay; no runtime is promised.
If too few qualify, use CLI `replay --top-k 10` or a larger pilot; never relax
scientific thresholds merely to fill the benchmark.

Freeze this checkout after preparation. Code changes invalidate execution against
the prepared inputs. Do not pull updates during a run; start a new output root
after changing code or scientific inputs. Existing completed outputs are not
overwritten. For a failed pilot, use a fresh output root or the CLI's identical
`evaluate --resume`; resumed execution is excluded from timing comparisons.

## One reference, six strategies, nested scales

After the pilot passes and resource usage is acceptable:

```bash
bash scripts/run_mdm2_matched_benchmark.sh full
bash scripts/run_mdm2_matched_benchmark.sh replay
```

`full` evaluates every held-out conformer in the common regular universe, with no
100,000-candidate early stop. It can be expensive: all templates, rejected schedules
and qualifying pose evidence are retained. It is explicitly launched, never run
by preparation. The fixed 100,000 cap applies to final candidate CSV ranking only.
Repeated `full` resumes an identical interrupted execution or verifies a completed
one. Complete reference evaluation is performed once for coverage, not once per
partition or scorer. Additional pilot and timing runs serve different purposes.

Replay derives 1/5/10/25/50/100% molecule-hash subsets from that reference, reranks
each subset and evaluates six trained orders plus 20 random orders per partition.
Conformers of each molecule stay together in scale selection. These are nested
subsets of frozen partitions, NOT independently rebuilt indexes at each size.
Their curves test routing scalability under fixed indexes; broader index-building
scalability requires separate rebuild and setup measurements.

| Output in `replay/` | Interpretation |
|---|---|
| `scale_summary.csv` | Actual molecule/conformer/qualified counts at each scale |
| `top_block_checkpoints.csv` | Discovery at Top1/5/10/20/50/all, including random trials |
| `routing_curves.csv` | Complete block-prefix discovery curves |
| `recall_at_work.csv` | Recovery within the same conformer-processing budget |
| `work_to_recall.csv` | Work needed for 90/95/99% reference recovery |
| `delivered_recall.csv` | Actual prefix-local reranking and output-cap recovery |

Report both conformer and unique-molecule recovery. Compare equal work, not only
equal block counts: block populations differ. Random trials must be summarized
with their variability, not just the best trial. Top-K values unsupported by a
small scale are omitted, not padded. Discovery means a reference candidate was
qualified in the visited blocks; delivery additionally requires that it survives
the prefix-local rank and export cap. These can differ.

## Real timing, not replay-time speedup

```bash
FRACTION=0.01 TOP_BLOCKS=20 REFERENCE_K=1000 \
  bash scripts/run_mdm2_matched_benchmark.sh timing
```

This deliberately executes 21 NEW searches: three complete subset searches and
three repeats of each of six Top20 methods, with rotated run order. Use a fresh
timing directory, same worker count and hardware, no competing tasks, and document
OS cache policy. Do not drop system caches automatically. Run a selected larger
`FRACTION` only after estimating cost; use the same fraction and output cap for
each paired comparison. `REFERENCE_K` must fit every full reference.

The summary accepts a speedup only when each repeat delivers at least 95% of the
specified reference Top-K AND the prefix output count matches the full output
count under the same 100,000 cap. Otherwise `accepted_speedup` is empty; the raw
time ratio is still diagnostic, not a successful matched-quality speedup. If Top20
fails, report that result and assess a larger fixed prefix separately. Choosing
that prefix after seeing this MDM2 reference is exploratory, not held-out tuning.

Timing includes input verification/scheduling, actual evaluation, ranking, CSV
writing and output sealing. It excludes MOL2 export, docking, historical partition
construction and priority fitting. Preparation duration is recorded separately;
upstream costs require their original receipts. Report cold-start total cost and
amortized repeated-query cost separately before claiming deployment advantage.
Replay-derived scales have no measured runtime unless explicitly evaluated.

The CLI also supports individual runs, arbitrary scale lists and explicit output
paths: `python scripts/benchmark_matched_blocks.py --help` and subcommand `--help`.
No experimental activity labels or independent target-level validation are added.

## Return for analysis

Return the small CSVs and `report.json` files from preparation, equivalence,
replay and timing summaries; include per-execution reports and protocols for timing
audit. Keep the large SQLite/chunk files on the workstation. First inspect recovery
at equal work, delivered recovery, and measured time at matched quality; do not
select a winner solely because Top20 returns many candidates.
