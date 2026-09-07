# E027 Dual-cocrystal overlap and docking-queue analysis protocol

Status: protocol locked before implementation

## Fixed input

- Accepted E026 aggregation directory containing two real WEE1 queries:
  8BJU/QT9 and 1X8B/824.
- E026 reported 12,834 union molecules, 955 shared molecules, 5,000 admitted
  molecules and 5,834 receptor-specific docking tasks.

## Hypothesis

The two queries are complementary rather than redundant. This predicts modest
Top-K overlap, substantial query-exclusive candidate sets, and a nonzero set of
shared admitted molecules that expand into two receptor-specific docking tasks.

## Analysis

1. Verify every aggregation output against its manifest SHA-256.
2. Recompute molecule counts per query and the all-query union/intersection.
3. Measure pairwise Top-100/500/1,000/5,000 intersection, Jaccard and overlap
   coefficient using within-query molecule ranks only.
4. For shared molecules, report Spearman correlation of within-query ranks.
5. Measure how often shape, unweighted color and anchored color choose the same
   conformer, separately for each query.
6. Recompute admitted support multiplicity, admission-reason counts and docking
   task counts per query/receptor.
7. Emit machine-readable JSON and a human-readable Markdown report.

## Confirmatory acceptance criteria

- Manifest hashes match all four input JSONL files.
- Recomputed union/intersection and admission/task counts reproduce E026.
- Every admitted molecule has exactly one task per supporting query/receptor.
- Pairwise Top-K counts are deterministic and intersections are nondecreasing.
- Raw cross-query Gaussian scores are never averaged or compared directly.

## Classification boundary

This analysis measures retrieval complementarity and queue construction. It
does not establish chemical diversity, binding, pose correctness, docking
accuracy or biological activity.
