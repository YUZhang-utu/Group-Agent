# Overnight handoff: 2026-09-21

User reports the workstation full-library validation is running and appears much
faster. Completion, actual throughput and final hit counts remain unverified.
Expected launch code: 1533883. Actual command, PID and output directory were not
provided. Do not start another job or update the running checkout.

## Active scope

- Library: 25,813,808 conformers / 8,318,351 source-grouped molecules.
- Query: 8BJU:QT9:A:601; ALL A0:HBD:30, A1:HBA:29, A2:HBA:39 and
  F12:pi_stacking:A:433:CG/CD1/CE1/CZ/CE2/CD2; minimum score 0.5.
- Joint eligibility: heavy_atom_ratio [0.7,1.3], maximum_extent_distance 0.35,
  minimum_feature_coverage 0.7. Exploratory rule, not biological acceptance.
- Suite intended settings: 22 NumPy workers, chunks of 2048, no Top-K/Top-N cap.
- Stages: initial controls, full-library rule evaluation, all matching molecule
  representatives, up to 128 final-hit reference checks. No automatic docking.

## Evidence already received

- Coarse audit: 10,032 checked, 288 retained, 9,744 rejected (97.129%);
  3.674 seconds, single CPU process, no seeds or Gaussian.
- Smaller hardware panel: 288 checked, 41 coarse survivors, all 41 rejected by
  seed feasibility. Reference hits zero. NumPy22 optimized median 1.079454 s;
  unfiltered NumPy22 median 1.575804 s. Positive retention not covered;
  scalability gate false. Different panels explain different rejection rates.
- Local overnight-suite regression: 331 passed, 3 skipped.

## Resume tomorrow

First inspect the existing run; do not rerun the launcher just to obtain status.
If the earlier automatic launch block created LATEST_RUN.txt, use:

```bash
BASE=/mnt/local/hand/yuzhang/aidd/e049-full-validation
OUTPUT="$(cat "$BASE/LATEST_RUN.txt")"
cat "$OUTPUT/suite-report.json"
tail -n 60 "$OUTPUT/run.log"
```

If launched directly with bash or with a different --output, use that actual
output directory instead. The pointer is not proof of which process is active.
A foreground launch may have no run.log; use the terminal output in that case.

Inspect full-library/progress.json while scanning. Once complete, inspect
full-library/report.json for full_coverage and per-stage conformer/molecule
counts, all-matches/representatives.jsonl for the complete deduplicated result,
and final-hit audit paths recorded in suite-report.json for reference agreement.
Check synthetic controls separately from real-library positive counts.

If incomplete but still running, leave it running. If interrupted, confirm the old
process is stopped before resuming the same output with unchanged code, selection,
workers and chunk size. Verified full-library chunks are reusable. Do not infer
absence of library hits from a zero-positive small panel.
