# Seed kernel and search-width audit

The current workstation job should finish on its original checkout. Do not pull
into that checkout, restart its workers, or resume its receipts with this version.
This change does not modify the running job. Reference generation and the 512 pair
cap remain the defaults. Batched generation is experimental and opt-in.

## Tomorrow: replay the same molecular panel

Use a separate checkout (or update the original checkout after its job finishes).
Activate the same workstation environment. RUN is the existing search directory
containing protocol.json, retrieval.json, selected-molecules.npy and
adopted-design/report.json. It is not the outer Chat project or page directory.

```bash
conda activate aidd-workstation
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export AIDD_ASSIGNMENT_BACKEND=numba
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export RUN=/absolute/path/to/execution/guided/guided/search
export AUDIT=/mnt/local/hand/yuzhang/aidd/seed-audit-$(date +%Y%m%d-%H%M%S)
python -u -m aidd_agent.seed_budget_audit \
  --run "$RUN" --output "$AUDIT" \
  --molecules 64 --workers 20 --chunk-conformers 32 --top 20
```

The 32-conformer audit chunks are deliberately small for load balancing on a
64-molecule panel; they do not replace the 4048-conformer production setting.
All stored conformers of the sampled retrieved molecules and all templates are
replayed. Sampling is by molecule, never by conformer. The original ranking and
MOL2 deliveries are not changed. No ANN search or full library pose rerun occurs;
collecting the sample's conformers does scan the catalog's molecule-ID arrays.

Eight variants: reference 512; batched 128/256/512/1024/2048; survivor 100/200,
each with a hard cap of 4096 pair seeds. PCA seeds are additional to the pair cap
and count toward survivor targets. Survival means passing the existing physical
mask; it does not imply distinct pose basins or biological correctness.

Read report.json for exact reference/batched-512 pose payload and ranking
equivalence; a failed equivalence makes the command exit with status 2. Compare
reference-top recall, molecule-template contact improvements, counts and timings.
The largest generated variant is the ranking reference, not ground truth.
Each variant/report.json includes every conformer-template's generated/surviving
counts and stop reason, including zero survivors. ranking.sqlite retains its
entire molecular ranking. Scores are contact-first, with Gaussian tie-breaking
and template RRF, as in the budget search.

The panel is selected from the ANN pool: no full-library recall claim is possible.
This first panel is exploratory. Repeat with a new output and seed, larger panels,
and independent chemistry strata before choosing production budgets. Sequential
variant timings include worker startup/cache effects. No automatic cap adoption.

## Kernel-only benchmark

```bash
python scripts/benchmark_seed_kernel.py --output "$AUDIT/kernel.json"
```

This alternates kernel order over 12 deterministic panels and checks exact seed
payloads. It excludes collision/scoring/I/O and is not an overall speed estimate.

## Explicit future budget-screen settings

The existing budget_screen CLI additionally accepts:

```text
--max-pair-seeds 1024 --seed-backend batched
--max-pair-seeds 4096 --survivor-target 100 --seed-backend batched --seed-batch 64
```

Use these only with the normal required batch/recommendation/output arguments
and a fresh output. Settings are sealed in protocol.json and cannot change during
resume. Chat defaults are unchanged; the audit command is the testing entry point.

## Block-tail input

Provide a CSV with block_id,molecule_id,score. Scores must come from a common
documented scoring and search budget; use contact score for contact-tail analysis.
Conformer duplicates within a block are collapsed to their best molecular score.

```bash
python -m aidd_agent.block_tail_audit --csv sampled-blocks.csv --top-k 5 > block-tails.json
```

The report includes max, top-k mean, effective k, unique sample size and optionally
high-score fraction when --threshold is explicitly supplied. Unequal samples are
flagged. This utility does not construct clusters, assign unsampled molecules a
score, prioritize a live job or reject blocks. Adaptive block retrieval and
pharmacophore/shape multi-route recall validation remain separate future work.
