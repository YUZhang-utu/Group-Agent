# E106: restore block-local consensus funnel and reduce repeated work

Protocol written before implementation validation, 2026-10-07.

Hypothesis: the reported block search cost is dominated by deterministic feature assignment; restoring the already adopted consensus threshold funnel reduces expensive pose work, while exact assignment acceleration and sealed block reuse reduce redundant computation.

User evidence: task 7fa0b402ba4b4a3a is running. Two blocks completed (37106 and 9792 conformers), third running. Bounded completed-chunk samples attribute 83.7% of recorded worker seconds to anchor assignment. This is not random sampling or a fresh benchmark.

New campaigns must use the adopted minimum_score, coarse_constraints, mandatory/alternative rules, geometry, pocket and minimum_pose_score through the same preselection compute path as consensus_funnel. Never invent conditions or treat optional weights as required anchors. Template passes are unioned by conformer; rank only actual survivors. Finish each started block before the 100000 survivor budget check; retain shortfalls.

Keep old requests on their recorded legacy mode. New request identity includes the new engine. Do not hot-update the active workstation or reuse old budget-mode chunks as threshold results. Shared reuse requires identical engine/code/query/library/member IDs/chunk schedule/backend, verified receipts and completed seals; reused elapsed time is not independent benchmark latency.

Validation: compare deterministic assignment arrays against the original solver including ties/zeros/rectangular matrices; compare threshold-mode explicit-ID output against the existing full-library compute path; test shortfalls, whole-block stopping, cache reuse/corruption/identity isolation. Run focused regressions and English guard. Timing fixtures are exploratory engineering measurements only; workstation speedup and current-target active retention remain unverified until separately measured.

## Local result
115 related tests passed, 1 skipped (optional dependency), English guard passed on 880 maintained files. Exact assignment checks passed including ties and zero weights. Reproducible synthetic solver measurements: about 0.97x for 8x12, 2.61x for 20x40 and 4.00x for 40x60; these do not establish end-to-end workstation speedup. See data/e106-assignment-check.json. Cleanup tests use isolated temporary files only. User authorized removing the old twelve-arm comparison; workstation cleanup/restart has not been executed from this Windows session.
