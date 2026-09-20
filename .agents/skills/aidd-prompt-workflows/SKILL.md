---
name: aidd-prompt-workflows
description: Configure and test the OpenAI-compatible AIDD planner and Project-scoped executor for protein data, AF3 and calibrated 3D search using structured English plans.
---

Read the [E038 setup guide](../../../to_human/E038_PROMPT_PROTEIN_AF3.md). Use
`scripts/run_e038_prompt_smoke.sh` for credential-free synthetic checks and
`scripts/run_e038_prompt.sh` for live planning/execution. A passing smoke test does
not establish live provider, protein API, AF3 or real-library success.

## Planning contract

The provider profile selects an actual endpoint/model and an API-key environment
variable. Never put the secret in configuration JSON, logs, plans or Git. JSON mode
can be disabled for compatible providers that lack it; local validation still
requires strict JSON and allowlisted actions with preceding typed dependencies.

Use `--plan-only` to inspect a plan. Normal execution may retrieve evidence and
prepare inputs; expensive AF3/search actions require the existing compute flag.
Apply authorization already provided in the session without redundant approval
questions. Genuine missing target identity or unsupported capabilities should yield
clear English clarification questions, not guessed tasks.

Keep summaries, clarifications, maintained prompts and documentation in English.
Do not modify biological identifiers/sequences to enforce prose language. Model
context contains the user prompt, capability schema and compact workflow skill
groups; it does not automatically ingest files or execute SKILL.md content.

## Project and recovery

Use the owned active Project and AI audit tables. The wrapper creates/reuses its
named Project; low-level `aidd_agent.prompt_workflow` accepts explicit context.
`ask` creates a new plan. Resume the printed existing plan with the `run` command;
never reconstruct a different plan into the same output directory.

Preserve plan seals, code/runtime fingerprints, stage receipts and reports. Report
`complete`, `blocked`, `failed` and `not_run` accurately. A failed external adapter
is not successful because its inputs were prepared. Return the relevant plan and
report locations, excluding credentials and sensitive user prompt content.
