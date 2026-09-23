---
name: aidd-prompt-workflows
description: Configure and test the OpenAI-compatible AIDD planner and Project-scoped executor for protein data, AF3 and calibrated 3D search using structured English plans.
---

Read the [E038 setup guide](../../../to_human/E038_PROMPT_PROTEIN_AF3.md). Use
`scripts/run_e038_prompt_smoke.sh` for credential-free synthetic checks and
`scripts/run_e038_prompt.sh` for live planning/execution. A passing smoke test does
not establish live provider, protein API, AF3 or real-library success.

## Planning contract

Use the [E043 guide](../../../to_human/E043_EXHAUSTIVE_CLASSIFIED_SEARCH.md) for
chat search, classified evidence, manual selection, crystal-pocket preparation and
Glide execution. Source-task actions are constructed by the trusted local chat
coordinator, not accepted as arbitrary provider-generated plans. Classification
requires no PLIP. Docking is a separate explicitly requested stage.

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

Select `--provider gpt` or `--provider deepseek` using separate environment keys.
See the [provider guide](../../../to_human/AF3_CONTAINER_AND_LLM_PROVIDERS.md).
Provider changes create fresh plans; they do not silently replace sealed plans.

## Persistent chat

For explicitly requested budget searches and molecule handoff, use the coordinator
`budget` and `budget_page` intents described in [E060](../../../to_human/E060_CHAT_BUDGET.md).
The user's combined search/export request authorizes both stages. Use contact-first
template ranks and preserve reserve pages. Do not route this request to the older
threshold selection/confirmation workflow. Full-library runtime and current-target
ANN recall remain measured outcomes, not guarantees from a configured budget.

Use [the conversational workbench](../../../to_human/E041_CHAT_WORKBENCH.md) when
users want ongoing dialogue, task progress or result lookup. The queue executes
existing Project-scoped plans; status/results are read from reports. Keep task IDs
session-scoped. Existing owned runs can be attached without new compute. A browser
refresh must not submit the scientific request again. Cancellation retains outputs;
resume remains bound to code/configuration and the original sealed plan.

The router receives bounded recent dialogue and task summaries; do not claim it
only receives the latest prompt or that it automatically sees all local files.

The E042 coordinator also sends compact anchor evidence and count summaries for
screening follow-ups. Use separate evidence, selection-preview and export turns;
the model must not invent thresholds or combine preview and export. Actual IDs,
counts and paths are supplied by local reports. See the
[screening handoff guide](../../../to_human/E042_SCREENING_TO_DOCKING_HANDOFF.md).
Do not present unavailable docking/new-target stages as completed capabilities.
