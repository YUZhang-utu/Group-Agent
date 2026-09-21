# Validate coarse survivors and the rule itself

## Overnight whole-library suite

The full suite is `scripts/run_e049_full_validation.sh --selection PATH --output DIR`.
It performs the 10000-row pre-run correctness audit, then processes EVERY library
ID under the same explicit joint rule with 22 NumPy CPU workers. No Top-K, Top-N,
max-chunks, or molecule export cap is supplied. Gaussian scoring applies to all
survivors of the declared coarse and seed-feasibility conditions, not every ID.

It validates full conformer/molecule coverage and verified chunk receipts, writes
the complete deduplicated representative list, and samples up to 128 final matching
molecules across that list for fresh reference reproduction. This final reference
audit is sampled; it is not a duplicate reference scoring of 25 million conformers.
Zero library hits remain explicitly unvalidated for positive retention.

Use `nohup bash scripts/run_e049_full_validation.sh ... > run.log 2>&1 &` in the
activated workstation environment. Inspect `suite-report.json`, `run.log`, and
`full-library/progress.json`. Final counts are in `full-library/report.json` and
the complete matching molecule list is `all-matches/representatives.jsonl`.
No SDF export or docking runs automatically. Duration is not guaranteed overnight.

For an interrupted run, stop any still-running original process, then rerun the
same command with the SAME output, code, inputs, worker and chunk settings.
Verified full-library chunks resume; small audits rerun into separate directories.
Do not update code while the job runs. A failed scientific check stops the suite;
zero-hit but equivalent exploratory panels do not prevent the requested full scan.

## Bounded survivor-only validation

Use the NEW selection preview report that contains `coarse_constraints`, not a
search, coarse audit or execution wrapper report. Its JSON kind must be
`selection_preview`. Reuse the same policy and unchanged source artifacts.

```bash
python -m aidd_agent.funnel_benchmark \
  --selection "$SELECTION_REPORT" \
  --output "$VALIDATION_OUTPUT" \
  --count 10000 --validate-survivors --no-gpu
```

This is a terminal entry point; no new chat command is claimed. It reconstructs
the deterministic 10000 spread IDs plus saved boundary/positive examples, saves
retained IDs, and evaluates ALL coarse survivors plus up to 64 rejected controls.
It does not run full-library refinement or impose a Top-K membership budget.
With the reported unchanged E048 inputs, the reconstructed panel should contain
10032 IDs and 288 survivors; verify rather than assume those counts.

The validation runs serial CPU batches to keep the experiment bounded and easy
to inspect. It is not a 22-worker throughput benchmark. The outer wall_seconds
describes the coarse phase only. Progress logs report the reference-positive count.

Inspect `report.json`:

- `validation.synthetic_control`: identity and Gaussian-selected query self-match,
  joint eligibility, anchor bounds and seed feasibility. These should pass for
  the current exploratory rule. Failure calls for inspecting the scoring/feature
  representation; do not silently weaken thresholds.
- `validation.reference_positive_count`: actual sampled library passes. Synthetic
  self-matches are deliberately excluded.
- `validation.checks`: all survivor arrays exact, checked IDs exact, zero false
  rejections. Any failed comparison makes the run failed.
- `validation.anchor_diagnostics`: individual and cumulative ALL passes across
  coarse survivors. Counts are pose-specific; do not sum them across conformers.

If no library positives remain, the run may complete with matching implementations,
but `positive_membership_validation` remains not covered. Inspect the per-anchor
breakdown, expand the exploratory panel if appropriate, and review poses before
changing the explicitly requested chemistry. No automated full-library launch.

The self-control uses the prepared query feature representation, not an independently
prepared SDF or experimental activity measurement. It cannot validate chemical
feature extraction, docking accuracy, or biological relevance by itself.
