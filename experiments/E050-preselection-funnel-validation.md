# E050: preselection funnel selectivity, retention and latency

Status: full-library entry implemented; workstation measurements pending.
Classification: exploratory. No workstation job is launched by this document.

## Question and hypothesis

Can the pre-human-selection 3D funnel substantially reduce expensive work and
end-to-end latency while preserving eligible candidates with at least one of the
specified crystal-derived anchor matches, together with their pose-specific
anchor combinations? The user selects important combinations after clustering.
Four-anchor ALL is not the early eligibility rule.

## Existing evidence and limits

User-supplied E049 report: 25,813,808 conformers / 8,318,351 molecules covered;
647,638 coarse survivors, all rejected by original-seed feasibility under ALL;
zero Gaussian-evaluated library conformers and zero final hits. Invocation time
895.4445416328963 seconds is not a corrected-ANY throughput benchmark.
Crystal-coordinate direction control and three synthetic rigid recovery cases
passed. Twelve selected conformers had no four-anchor passing original seeds or
expanded anchor-directed seeds. These results do not measure ANY selectivity,
whole-library absence of feasible poses, or biological performance.

## Locked pilot semantics

- Use the selected QT9 query and four existing anchor IDs; at least one passes.
- Initially keep score threshold 0.5 and the recorded whole-ligand eligibility
  thresholds for comparability. Treat those coarse thresholds as provisional
  eligibility policy, not mathematically safe bounds on anchor membership.
- Keep original seed parameters fixed for the primary acceleration comparison.
  Expanded search is a separate sensitivity experiment, not the reference for
  claiming unchanged original-seed acceleration.
- Reference annotates all declared original seeds for pilot conformers without
  early anchor pruning. Preserve per-pose anchor bitmasks and scores, not only
  the Gaussian winner. Molecule unions must not imply simultaneous contacts.
- Evaluate coarse-policy exclusions separately: score stratified excluded rows
  to quantify loss of anchor-positive candidates outside the eligibility policy.
- Clustering retains membership and links to pose evidence. Representative-only
  display is not deletion of members or an improvement in chemical selectivity.

## User-authorized direct full-library sequence

The user explicitly rejected pilot prerequisites and requested the complete
library now. This supersedes the earlier small-to-large scheduling proposal.
Local implementation tests are not workstation scientific pilot runs.

1. Fixture-test ANY pruning, alternative-pose retention, molecular aggregation,
   full coverage, sealed-chunk resume and stage counters.
2. Run preselection_full on every catalog ID, 22 workers, 2048-row chunks, with
   no max-chunks, Top-K, Top-N or pre-run sample. Use a fresh E050 output.
3. Keep original seed parameters. Compute Gaussian conformer winner scores for
   ranking, and exact anchor assignments for all optimistic ANY-passing seeds.
   Keep a real representative for every observed exact anchor bitmask per molecule.
4. Group all retained molecules by exact non-stereochemical Murcko scaffold;
   keep all members and combinations. Acyclic molecules use full connectivity;
   reconstruction failures remain explicit singleton groups. This is structural
   grouping, not fingerprint-distance or 3D similarity clustering.
5. Report full coverage, measured stage counts and timings and human-selection
   tables. No unpruned whole-library equivalence or speedup claim is made by this
   single run. Candidate recall outside coarse eligibility remains unmeasured.

## Required measurements

For each size/type-coverage/extent/anchor-bound/seed-feasibility/scoring stage:
input and output conformers, unique molecules, stage and cumulative rejection
fractions, wall time, seed count, exact score evaluations, and reference positives
lost. Attribute sequential rejection reasons in the actual execution order.
Also report counts by exact pose anchor bitmask and by number of matched anchors,
deduplication time, clustering time, cluster/member counts, output serialization,
peak memory, total elapsed time, and matched-concurrency speedup.

## Acceptance and interpretation

No fabricated per-stage rejection target. A stage is useful when its cost is
justified by measured downstream savings without unreported candidate loss.
Require zero observed false rejections and preserved requested pose evidence
against the declared pilot reference, with real-library positives present.
No-positive agreement cannot validate positive retention. Coarse eligibility
losses remain a separate policy tradeoff; do not exclude them from reporting by
definition. Full preselection acceptance requires clustering and human-selectable
anchor-tagged members, not merely an ANY scoring benchmark. Projected full-library
time is an estimate until measured. No docking or affinity validation is claimed.
