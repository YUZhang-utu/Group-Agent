# E109: matched evaluator, routing replay and measured prefix execution

Protocol before implementation/fixture tests, 2026-10-14. E108 outcomes have already
been inspected; this is exploratory method development, not preregistered evidence.

Hypothesis: E094/E095/E096 and ChemPLP/EquiScore priorities reduce the work and time
needed to recover a fixed global reference Top-K under one frozen evaluator.
Larger libraries may increase the advantage; this is tested, not assumed.

Freeze the adopted eleven-template query (shape/color, contact and 5C5A pocket),
current threshold branch, NumPy assignment backend, pose generation and positive-
contact/per-template contact-Gaussian/RRF(k=60) ranking. Keep existing 7NA2_A block
priorities as an explicit routing hypothesis. No new docking or priority training.
Do not change weights, adopt another receptor or enable new constraints.

Prepare one source-verified common regular population, all three memberships,
canonical artifact/molecule IDs and all six full block-priority lists. Exclude all
molecules present in the priority-training sample by default, across conformers.
Record exact counts and the excluded special/chemistry-review scope. Deterministic
molecule-hash fractions define nested scales; report actual counts, not requested
nominal sizes. Query-level holdout remains absent for the current MDM2 experiment.

Evaluate one complete chosen universe through the existing threshold evaluator.
Preserve scheduled rejected rows as well as scored poses. Do not reuse prior E106
result caches. Derive smaller-scale reference ranks from the saved per-template
poses, reranking within each scale without rescoring. No historical ANN output
serves as exact truth. Existing scientific files remain read-only.

Replay trained and random block orders. Report cumulative discovery at Top1/5/10/
20/50/all, fixed work fractions and work to 90/95/99% recall, for conformers and
canonical molecules. Separate routing discovery from final output: optionally
recompute local prefix RRF and report delivered reference recall with the same
output cap. Global-reference rank must never be used as the deployed prefix rank.

For actual timing, execute fresh full or chosen block-prefix ID sets with the same
workers, chunk size, evaluator, query and aggregation. No candidate-count early
stop is imposed in this benchmark. Prefixes are fixed before execution, not chosen
by peeking at reference outcomes during deployment. Report protocol/ID hashes,
resume/cache status and timing phases. Pair at least three trials before reporting
time ratios. Require a successful complete-all-block equivalence check on a bounded
pilot; all-block and full candidates must match after deterministic aggregation.
Replay wall time is never search speed. Derived-scale references have no measured
baseline time. Setup, priority scoring, storage and amortization are separate.

Validate identity mapping, partition conservation, train/test exclusion, nested
scales, ranking equivalence, duplicate/missing priorities, deterministic random
controls, full/prefix containment, resume sealing and ineligible timing exclusion.
Only synthetic tests execute locally; workstation science requires explicit CLI.

Implementation validation, 2026-10-14: 13 E109 synthetic tests and 18 related
historical-routing/configuration tests passed (31 total). The English-content
guard passed. Scientific engine calls are mocked in E109 fixtures; the real
workstation full/all-block pilot remains mandatory. No new MDM2 scientific result,
activity enrichment, winning strategy or measured speedup is claimed.
