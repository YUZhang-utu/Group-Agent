# E127: reliable planning, task receipts and failure recovery

Date: 2026-11-10. Supersedes the E126 task-failure behavior.

## Confirmed incident

The user supplied task `0a928af2a73a479b` and its status response:
`failed / Invalid plan fields / No completed scientific report is available yet`.
The error is raised by the workflow planner's top-level schema check before a plan
is saved or a scientific engine starts. This is not evidence of a bad CIP2A sequence,
invalid SMILES, GPU exhaustion or failed AF3 inference. The original invalid model
response was not supplied, so its exact missing/extra fields cannot be reconstructed.

The conversational trace separately showed a successful queue submission followed
by two malformed final-answer attempts. Those were a second failure surface; they
did not explain the background plan failure.

## Changes

- State the exact workflow-plan contract in the schema and planner instructions:
  `version`, `summary`, `clarifications`, `steps`; each step has `id`, `action`,
  `params`. Report missing/extra/allowed fields instead of a generic field error.
- Retry invalid model JSON/schema/provenance at most twice after the first attempt,
  using validation feedback and the same original request. No scientific step or
  executable run is created before validation. Authentication, connection and
  configuration errors are not blindly retried. Repairs do not retry inference.
- Return a deterministic queue receipt after a workflow submission. Do not invoke
  the conversational model again merely to describe an already submitted task.
  Queue admission does not claim that an input file exists or inference has started.
- Create task logs before planning. Record phase-specific failure diagnostics and
  display the reason directly in Chat. Task-report lookup returns stored errors
  and paths even if no scientific report exists. Preserve credential redaction.
- Surface real step errors after worker execution, including configuration blocks.
  Show current steps, prepared input paths and final model paths in task cards.
  Long request text is expandable; completed tasks offer results rather than resume.
- Return discoverable AF3 input/provenance and prediction-manifest paths. These
  references are also visible through `/results` and the owned artifact tools.
- Preserve complete unambiguous labelled protein sequences and explicitly labelled
  seed lists. Reject model-generated sequence truncation, prohibited constructs or
  altered seed lists. Mark the current user turn so later seed changes can override
  an earlier request. Existing generic/verified-protein modes remain available.
- Label the isolated `qa_chat_features.py` page as a synthetic interface test.
  It must not be mistaken for the configured Agent that runs scientific tasks.

The repaired planner is shared by the existing scientific task workflows. This
does not introduce new engines or complete the previously pending generic MOL2
export adapter/shared-MSA batch Chat dispatcher. Runtime installation and scientific
quality checks remain separate from successful task orchestration.

## Recover the existing CIP2A task

In the existing repository and AIDD environment:

```bash
git pull --ff-only origin feature/structure-guided-chat
python -m pip install -e . --no-deps
```

Restart the configured Agent server using the same storage, runtime and provider
arguments, including `--allow-compute` for actual AF3 execution. Restart when no
other scientific task is running. Open the newly printed access URL, return to
the original conversation, then enter:

```text
/resume 0a928af2a73a479b
```

This existing task failed before a plan was created, so resume retries planning
under the new code with the recorded request and original task ID. Do not create
another copy of the prediction request. This instruction is specific to the
reported pre-plan failure: partially executed sealed workflows may reject resume
after a code change and require separate inspection.

Use `/status 0a928af2a73a479b` or **Show results** to inspect it. A completed input
step should expose `input_json`; the prediction step can still report a specific
runtime/compute failure. A queued receipt is never a completed AF3 prediction.

## Verification

The test fixture preserves the supplied 905-residue sequence, exact ligand SMILES
and seeds `[47, 48, 49, 50, 51]`. Sequence SHA256:
`30683fdc89a524c0b0a383343a3e17ba2c54b4dfe8d5d37b22856a175d97e584`.

- 119 tests passed across planning, Chat, domain tools, streaming/uploads, generic
  query/AF3 and structure-chat suites in 33.68 seconds.
- 67 additional container, prediction, guided/budget, chain, import and AF3-batch
  cases passed in 6.46 seconds.
- After adding original-task recovery, the 27-case reliability/HTTP suite passed
  in 8.24 seconds. This overlaps earlier suites and adds one new test: 187 distinct
  Python cases passed across the targeted runs.
- Four Markdown checks and JavaScript syntax checks passed. Browser visual
  acceptance was not performed in this update.

End-to-end local checks include the actual Chat queue and worker subprocess,
malformed-plan repair, strict ligand parsing, exact AF3 JSON, retained input
artifacts, and resuming a failed planning job without changing its ID. GPU execution
is intentionally blocked in these fixtures; no live DeepSeek/GPT call, production
AF3 inference, docking or library scan was performed. The workstation must confirm
the next runtime stage; these checks do not guarantee every provider response or
every scientific run will succeed.
