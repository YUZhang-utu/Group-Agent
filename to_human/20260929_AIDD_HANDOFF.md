# AIDD handoff: 2026-09-29

## User instruction and immediate priority

The user reports that the workstation macrocycle routing task is still running.
Wait for it to finish before the user pulls the new code. Do not stop, restart,
replace, or rerun that job, and do not request a code update during its execution.
No workstation process has been inspected directly from this local workspace.
Resume by reviewing the finished job's receipt, not by launching another build.

## Workstation job: E091 calibrated boundary routing

Module: aidd_agent.calibrated_boundary_routing.
Last published implementation before E093: commit 5b9f25c.
Guide: to_human/E091_CALIBRATED_ROUTING.md.
Expected paths based on the provided run instructions (confirm from the receipt):

- Build: /mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone
- Candidates: same parent, backbone-boundary-candidates-v1
- New routing output: same parent, backbone-calibrated-routing-v1

The user confirmed a special-set engineering target of at most 1% of the
24,663,736 admitted conformers: at most 246,637 special conformers. The separate
1,150,072 source chemistry-review records are outside this denominator and must
not be forgotten or represented as resolved. Current full-run outcome is unknown.
The <=45 degree gate has 1,262,619 angular candidates, leaving a minimum121,720
special conformers before property/steric/backbone compatibility checks.
At least1,137,702 boundary assignments are required to reach the1% goal.
Do not loosen chemistry or chirality to force the target.

E091 uses at most64 distinct-molecule references per target, at least32 required,
50/25/25 fit/calibration/check splits, and at most12 prototypes. Membership is a
provisional sidecar; ordinary cis/trans classes and old databases remain intact.
It is slow because it verifies and profiles the full candidate population,
not only the300-member E089 review panel. No automatic interrupted-run resume.

After completion inspect report.json for status, special_fraction,
engineering_target_met, additional_assignments_needed_for_target, special reasons,
model availability/check failures, source provenance and stage timings. Inspect
boundary_assignments.sqlite and routing_models.json if reasons need diagnosis.
Verify conservation of ordinary, assigned-boundary and special populations.
If the target is missed, identify the actual limiting channel or reference
shortfall before proposing modifications. No production block rejection, docking
or screening-recall validation has been established.

## Completed locally and published: E093 receptor advice

GitHub: https://github.com/YUZhang-utu/Group-Agent
Branch: feature/structure-guided-chat
Implementation commit: 9515c31 (successfully pushed).
Guide: to_human/E093_RECEPTOR_ADVICE.md.
Validation receipt: to_human/E093_VALIDATION.json.
Protocol: experiments/E093-receptor-advice.md.

The existing pockets step now delivers concise receptor-selection advice:
qualified chains, independent PDB entries, scaffold-series proxies, provisional
single receptor for insufficient/unresolved evidence, otherwise supported partition
representatives or a coverage set with recommended count and coverage tradeoffs.
The user chooses the receptors. Detailed matrices remain supporting artifacts.
Weighted normalized shape fusion is versioned; old sealed reports are unchanged.
Individual receptor_options, adoption, viewer selection and local consensus are
connected. Coverage neighborhoods are not physical states or global mandatory
interaction intersections. Missing ligand preparation is explicitly marked.
Domain-agent instructions use the existing literature_search tool for cited
interpretation; live LLM/literature behavior has not been accepted on workstation.

Local checks:87 related tests passed; final budget fallback13 focused tests passed;
English-content guard passed. Exploratory MDM2 replay:67 chains/59 PDBs/25 generic
scaffold proxies, no robust partition resolved. Under current provisional scales,
8 references cover72.3%,21 cover95.2%;5 chains remain uncovered and retained.
This is not21 biological states or a validated requirement for docking21 receptors.
Several representatives still require ligand preparation before consensus.
Local review: data/e093-mdm2-final-advice/report.html.

## Next session order

1. Read this file, research-state.yaml and the latest research_log.md entries.
2. Ask for/read the completed E091 report, or report the job as still running.
3. Analyze the actual special-set outcome and integrity before changing routing.
4. Only after that job ends, the user plans to pull feature/structure-guided-chat.
5. Restart Chat when ready to use E093; perform workstation/live-provider acceptance.
6. Keep PLANTS/MDM2 and the user's N-E workflow as later evaluation steps;
   affinity model remains undecided. Do not imply they have executed.

Suggested continuity prompt:
Read to_human/20260929_AIDD_HANDOFF.md and research-state.yaml. Continue from the
E091 workstation routing result; preserve the running job if it has not finished.
Do not rebuild the library or infer that the1% special-set goal has been achieved.
