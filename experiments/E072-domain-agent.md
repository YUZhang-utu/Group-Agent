# E072: grounded multi-step domain agent

Implementation and acceptance protocol, 2026-09-25.

Hypothesis: one bounded observe/tool/result loop over existing validated adapters
supports compound research requests better than one-shot intent routing. This is
an engineering hypothesis, not evidence of superiority to Codex.

Implement session/task discovery, owned report reads, configured-library lookup,
literature search with provenance, existing workflow execution and live viewer
feedback. Preserve compute flags and queue receipts. Persist every decision and
tool receipt, reject duplicate mutations, and never equate queued with completed.

Tests: compound report/library/literature synthesis; tool failure recovery;
ownership/path escape; source citation validation; duplicate task prevention;
pending viewer/job behavior; restart-safe audit; old slash-command regressions.
Use deterministic planner fixtures to verify orchestration separately from live
provider acceptance. Document untested workstation/desktop boundaries explicitly.

Live benchmark: paraphrases of at least ten real domain tasks on the workstation,
with completion, correctness, evidence traceability, time, tool count and unwanted
job launches. Compare the old router and new agent on the same task/data snapshots.
Do not claim a benchmark win without these measurements.
