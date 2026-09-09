# E031 General fast chemical reranking and docking-handoff protocol

Status: protocol locked before implementation

## Clarification before implementation: final retrieval output

The final pre-docking layer answers one query-centric question: how well does
an already retrieved and rigidly aligned candidate reproduce the key
interactions observed for the co-crystal ligand? It does not search new poses.

For every query interaction anchor `i` and compatible candidate feature `j`,
evaluate a bounded `[0,1]` match from feature type, spatial agreement, and
direction/axial agreement. A deterministic one-to-one assignment prevents one
candidate feature from satisfying multiple distinct query anchors:

```text
interaction_match_score = sum_i(weight_i * best_assigned_match_i)
                          / sum_i(weight_i)
```

Persist each per-anchor contribution, assigned candidate feature, and unmatched
anchor. Extra candidate features do not erase a reproduced query interaction.
Shape and ordinary color remain earlier evidence lanes and are not mixed into
this score with arbitrary uncalibrated weights. A later combined ranking would
require held-out calibration.

The interaction schema is target-independent. Direct hydrogen bonds are the
first implemented class; salt bridges, aromatic stacking, cation-pi,
hydrophobic contacts, and metal coordination require explicit typed extractors
and validation rather than WEE1-specific rules.

## Objective

Preserve a general, fast, effective 3D retrieval path across targets and
chemotypes. WEE1 is the first validation case, not an algorithm-specific rule.
Receptor-dependent conformational relaxation belongs after retrieval in the
docking/pose-preparation stage.

## General stage boundary

1. One-time offline companion construction is amortized across all queries and
   excluded from warm query latency.
2. FAISS and reusable pharmacophore postings define broad recall.
3. Existing all-candidate PCA Gaussian and fixed Top-N rigid Gaussian remain
   the only pose-search operations in pre-docking retrieval.
4. Key-interaction coverage is evaluated once under each already selected rigid
   objective pose. It adds no new rigid or torsion seeds and emits one
   explainable score plus per-anchor assignments.
5. Receptor vdW is evaluated once only for the small, caller-capped docking
   handoff/export set. It is an annotation and cannot alter retrieval admission.
6. Torsion beams, rigid micro-seeds, local minimization, and docking are
   downstream operations whose latency is reported separately.

## Hypothesis

Mmap-backed candidate directions can add chemically relevant evidence to the
existing fixed rigid-refinement set with at most 10% incremental refine latency,
because they reuse candidate IDs, feature coordinates, and transforms without
expanding the pose search. Separating receptor-dependent relaxation preserves
front-end speed and makes the architecture general across targets.

## Performance gates

- Disabled chemistry path remains byte-identical to the accepted baseline.
- Enabled directional annotation adds no candidates, seeds, or transforms.
- Incremental directional wall time target: <=10% of the accepted rigid-refine
  wall time on identical candidate IDs and hardware; report absolute time too.
- No full-library receptor-distance calculation.
- Docking-handoff vdW work is linear only in an explicit small task cap and is
  excluded from retrieval latency.
- Peak RSS and mmap reads are measured separately for cold and warm execution.

## Effectiveness gates

- Preserve baseline protected recall and every original objective/rank.
- Store directional overlap as an additional named, reversible evidence lane.
- Measure rank correlation and Top-K overlap before enabling any reranking.
- Calibrate usefulness on multiple queries/targets and held-out enrichment or
  redocking benchmarks; two WEE1 structures cannot establish general benefit.
- Drop or leave direction annotation disabled if it fails to improve held-out
  early enrichment robustly.

## Implementation sequence

1. Batch-read companion feature memberships/directions only for rigid-refined IDs.
2. Apply each stored rigid rotation to directions without translation.
3. Build a small anchor-by-feature compatibility matrix and perform
   deterministic one-to-one assignment once per retained objective pose.
4. Persist the normalized score, assignments, unmatched anchors, and timing
   without changing baseline result arrays.
5. Benchmark identical IDs with chemistry disabled/enabled.
6. Separately expose single-pose vdW annotation in the docking handoff.

## Prohibited in the search lane

- No 13-seed micro-relaxation.
- No torsion beam.
- No receptor-specific hard clash filter.
- No WEE1-specific threshold or conditional branch.
- No uncalibrated weighted sum replacing accepted Gaussian objectives.
