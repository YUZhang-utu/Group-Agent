# Target project direction, knowledge and execution contract

Status: proposed architecture and operating specification, 2026-11-11.
This document does not mean the proposed dossier, context builder or execution
gates are already implemented. E131 remains the current runtime baseline.

## Primary purpose: keep work aligned with the user's intended outcome

The user's clarification makes direction control the primary purpose of the
dossier. Shared knowledge, context compression and task history support that
purpose. A technically successful task can still be the wrong task.

Pin an original mission at the start: the user's original request/source message,
intended outcome, deliverables, acceptance criteria, scope, exclusions and fixed
scientific requirements. Separate that mission from the agent's chosen methods.
Do not continually rewrite the mission to match whatever work has been performed.
Preserve the original baseline and a current authorized revision. Explicit later
user changes update the current revision with provenance; the baseline must not
prevent the user from changing direction.

Every task must name a mission requirement ID and explain its concrete contribution
to that requirement. Necessary preparation and investigation of blockers may be
dependencies, but must retain that link and a stopping condition. New interesting
findings belong in a backlog until they serve the active mission or the user
authorizes their pursuit. Finishing more tool calls is not progress unless they
satisfy a deliverable, resolve a needed dependency or close a relevant uncertainty.

At task planning/dispatch, after results, on resumption and before selecting the
next task, check alignment against the pinned authorized mission, not just against
the previous task or latest generated summary:

1. Which original/current authorized requirement does this action serve?
2. Are target, construct, assembly, chemistry, library scope and requested outputs
   still consistent with that requirement?
3. Does the selected method produce the needed evidence/deliverable, and are we
   reusing appropriate completed work instead of pursuing an unrelated branch?
4. Have an error or missing dependency led us to weaken acceptance criteria or
   substitute an easier but different task?
5. Is this a method adjustment within existing authorization, or a mission change?

Routine implementation choices and within-scope repairs proceed autonomously.
Use authorization already present in the conversation; do not repeatedly ask the
user to approve the same scope. A change to goals, fixed scientific requirements,
deliverables or acceptance criteria requires user direction if not already
authorized. Record a concrete proposed change and its reason; stop only dependent
work while continuing independent authorized work. Do not let the agent approve
its own proposal or treat a generated interpretation as a user instruction.

Examples for this AIDD project:

| Intended work | Direction drift to catch | Correct handling |
| --- | --- | --- |
| Target -> experimental ligand reference -> local-library candidates -> MOL2 | Ending after a literature review or reference download | Record completed prerequisite; retain screening/export as outstanding deliverables |
| Score every configured library conformer | Silently substituting a retrieved subset after a slow run | Preserve the coverage requirement; report the blocker or propose an explicit scope change |
| AF3 dimer with supplied sequence and exact seeds | Running one chain or default seeds because configuration is easier | Reject the mismatched plan; preserve inputs and fix or report the concrete blocker |
| Re-export a saved candidate set | Launching another search, docking or model-training project | Reuse the saved ranking/poses; keep unrelated ideas outside the active plan |
| Produce a scientific prediction | Declaring an interface fixture or queued receipt sufficient | Keep the scientific deliverable incomplete until actual required evidence exists |

Record completion by requirement: satisfied with evidence, outstanding, blocked,
or explicitly removed by the user. A blocked requirement is not silently removed.
The mission summary and requirement status should remain small and always present
in task context; detailed history and domain evidence can be retrieved as needed.

### A long-running project needs requirements, a plan and an execution ledger

Treat this like a software project whose implementation must continue to satisfy
its requirements over many work sessions. The dossier has three linked layers:
the user-authorized requirement baseline, the dependency/milestone plan, and the
actual execution ledger. Requirements describe the desired outcome and limits;
the plan describes the current method; the ledger records what really happened.
Do not confuse a planned step with executed work or an executed step with an
accepted deliverable. On a new chat or restart, reconstruct the next ready work
item from these records instead of replanning the whole target project from scratch.

### Check reuse before creating another scientific task

Each proposed calculation needs a stable work signature based on operation,
project/branch scope, exact input content hashes, scientific parameters (including
seeds), query/reference identity, library version, relevant tool/runtime version
and requested output semantics. Wording changes in a prompt or a new chat must
not create a different scientific identity by themselves. Use stage-specific
signatures: changing an export quota does not change an already computed full
ranking. Required reproducibility checks must still hold when outputs are reused.

| Prior work | Next action |
| --- | --- |
| Equivalent task queued/running | Attach to its progress; do not submit another copy |
| Equivalent task complete, outputs present and valid | Reuse its artifacts and link the source task |
| Required ranking exists; only selection/export differs | Create the downstream selection/export task only |
| Task failed/interrupted and checkpoint is compatible | Resume the existing task under the established recovery rules |
| Submission outcome is uncertain | Reconcile the receipt before any resubmission |
| Inputs/parameters/library version changed | Record the delta; rerun only affected stages and dependent descendants |
| Prior outputs are missing or invalid | Record the failure of reuse; recover outputs or justify recomputation |
| User requests an independent replicate | Record that intent and replicate identity; do not deduplicate it away |

Record the reuse decision and consulted task IDs before dispatch. A new task must
state why an existing result is insufficient. Scientific freshness and parameter
compatibility matter, not just matching filenames or finding a complete flag.
Reuse links must preserve provenance and never copy a prior result into a new
task as if a new calculation had occurred. Existing controls already protect some
submission paths; a unified semantic reuse guard across all operations is proposed,
not implemented by this document.

## 1. What currently shares information

Each conversation has an immutable project_id. Task lookup, uploads, imported
results, candidate sets and notes use that project boundary. The agent can read
task reports and search earlier project conversations without rerunning work.

The current domain loop supplies the latest 10 conversation messages, each cut
to 1500 characters, plus project memory. Project memory includes the latest 20
tasks (request excerpts limited to 350 characters), latest 20 notes, workflow
summaries and candidate-set summaries. Workflow and candidate-set lists do not
yet have a shared total context budget. The loop also includes up to six recent
tool-evidence items, individually compacted, and a tool schema. Its 200000-character
request ceiling is an input limit, not a validated context-quality target.

These mechanisms preserve access and provenance, but they do not establish a
curated knowledge base. An old constraint can disappear from the recent-note
window. A recent statement may be less relevant than an old verified artifact.
An execution marked complete does not establish biological validity. Free-text
notes are recorded statements, even when their kind is finding.

## 2. One readable dossier, backed by structured records

Give each target project a visible `PROJECT.md` using
[the dossier template](templates/TARGET_PROJECT_DOSSIER.md). It is a concise
working record, not a cumulative chat transcript. It answers:

First: are we still pursuing the originally intended outcome, as explicitly
revised by the user, and how does the next action advance it?

1. What exact target, construct and assembly are we studying?
2. What is the user's objective and which constraints remain active?
3. Which evidence and assets have been checked, and at what level?
4. Which task branch are we following, what has happened, and what is next?
5. What remains unknown or conflicts with earlier information?

The proposed authoritative state is structured project records, extending the
existing database with versioned constraints, assets, decisions and work items.
`PROJECT.md` is a deterministic rendered view of those records, not a second
independently edited source of truth. A user edit should be imported as a proposed
patch, validated and revisioned before replacing canonical state. Until that
feature exists, the template is a manual planning artifact only.

Keep full reports, structures, rankings and exact scientific inputs in their
existing local artifacts. Store durable artifact IDs, paths, hashes, task lineage
and verification dates in the project records. Keep historical revisions and
events in an append-only ledger, outside the compact current dossier. Both views
should be reconstructable without a surviving model conversation.

## 3. Separate claims, decisions and evidence

Every knowledge entry needs an ID, type, scope, source, revision and status.
Scope includes the project and, where relevant, construct, assembly, binding site,
reference ligand, library version and workflow branch.

| Entry type | Required provenance | What it authorizes or establishes |
| --- | --- | --- |
| User constraint | Original message ID and explicit decision | Task requirements; not a scientific measurement |
| Input identity | Immutable input artifact, hash and origin | Exact supplied or retrieved material; distinguish user-supplied from database-verified |
| Observed result | Task/report ID, field pointer, units and output hash | The recorded measurement or runtime observation only |
| Interpretation | Linked observations and stated uncertainty | A proposed explanation; not an independently established fact |
| Hypothesis | Rationale and a proposed test | Work to investigate, not a finding |
| Decision | Alternatives, selected option, author and rationale | Which branch to follow and under what scope |
| Open issue | Missing/conflicting evidence and affected tasks | An explicit dependency that may block dispatch |

Execution status and scientific validation are separate fields. A queued AF3 task
is not a prediction. A completed prediction is not a validated binding pose.
A report imported without payload checks is report-only evidence. A similarity
ranking is not an affinity ranking. An attached task ID alone does not validate
the prose of a note: the referenced field must support the associated claim.

Do not use a universal "latest statement wins" rule. A new user decision may
supersede an older constraint, with an explicit revision link. It cannot rewrite
a measured result. Conflicting evidence remains visible until resolved. A decision
for one construct/site must not silently govern another branch of the same target.

## 4. Assemble a small task-specific context

Build a fresh context for each task rather than carrying a growing transcript:

1. Current user request and exact requested operation.
2. Compact project identity, active constraints and unresolved blocking issues.
3. The selected task card: objective, branch, inputs, parameters and dependencies.
4. Relevant evidence excerpts with source IDs and verification status.
5. Relevant decisions and a small conversational tail if it resolves intent.
6. Only the tool contracts needed for the selected operation.

Always retain applicable constraints by their active status, not by recency.
Always retain the authorized mission revision, its acceptance criteria and the
current task's requirement linkage. These take priority over conversational history.
Retrieve evidence first by explicit task/asset/branch references; then use metadata
or text search if necessary. A semantic retrieval component is optional and must
not replace deterministic ownership, identity and version checks. Retrieval
relevance is not evidence verification.

Initial engineering budget proposal: reserve approximately 1500 tokens for project
identity/constraints, 2000 for the current task/decisions, 4000 for relevant evidence
and 1000 for conversational context. These are adjustable design targets, not
benchmarked optimal values. Count the provider's complete request, including tool
schemas and current input, with response headroom. Record included/omitted IDs and
the reason for inclusion. Never silently discard a mandatory constraint to fit.

Long protein sequences, SMILES tables, PAE arrays and full rankings should be passed
to executors as validated typed artifacts, with hashes, rather than repeatedly
summarized or copied by the LLM. Do not truncate scientific inputs to save tokens.
If essential material does not fit, narrow the task or retrieve a precise field;
report what is unavailable. A missing value stays unknown, never a guessed default.

## 5. Follow a target workflow with explicit branch choices

These stages define work, not a promise that every target has experimental
complexes or every downstream adapter is configured. Skip only with a recorded
reason; reuse valid artifacts when the inputs and requested purpose match.

| Stage | Work | Required exit evidence |
| --- | --- | --- |
| Scope | Record objective, target label/accession, species, construct and monomer/dimer composition | Identity provenance, exact input asset and unresolved conflicts |
| Inventory | Inspect prior tasks, local assets, runtime/library configuration and prior failures | Reuse candidates and missing dependencies; no new calculation implied |
| Route | Select experimental-reference screening, supplied-ligand screening, AF3 or another supported branch | Chosen route and requested deliverable |
| Experimental reference branch | Retrieve complexes; assess structure, ligand identity, sites and diversity | Exact PDB/assembly/chains/ligand instance; resolve distinct sites or references |
| Supplied-ligand branch | Inspect attachment/SMILES; preserve chemistry; prepare missing 3D conformers | Input/output hashes, actual conformer/query IDs and rejected-record report |
| AF3 branch | Prepare exact sequences, chain copies, ligand chemistry and seeds; check configured runtime | Saved input first; then real task status, model paths and recorded confidence |
| Library search | Select the actual query, library version, scoring route and output quota | Coverage counters, complete saved ranking and explicit scientific limits |
| Candidate management | Save/reuse ranking; apply requested filters, labels, set operations or scaffold grouping | Named set revision, parent sets and retained molecule IDs/counts |
| Delivery | Export original/posed MOL2 and names; verify identity/provenance | Actual output paths, hashes, counts, shortfall and failures |
| Next experiment | Review uncertainty and propose docking, AF3 comparison or another supported task | A new explicit task card with acceptance criteria; no affinity claim from similarity |

AF3 is a parallel/optional route, not a mandatory precursor to a cocrystal-based
search. The full-catalog Gaussian route and ANN-based budget route remain distinct.
Finding several valid ligand references does not authorize choosing the first.
Re-exporting or filtering a saved candidate set does not authorize another scan.

## 6. Task protocol and state updates

Before dispatch, persist a task card with:

- Goal and requested outputs; project, branch and parent task IDs.
- Mission revision, requirement IDs, expected contribution and alignment verdict.
- Exact input artifact IDs/hashes and selected construct/site/reference versions.
- Parameters, configured tool/runtime, compute authorization and output location.
- Execution completion criteria, scientific QC criteria and stop conditions.
- Project revision used to build the plan and any unresolved dependencies.

Dispatch only after deterministic checks verify required fields, ownership,
hashes, supported operations and applicable constraints. If the project revision
changes before dispatch, re-evaluate affected requirements. Running work keeps its
original snapshot; a later decision does not retroactively rewrite that task.

On a task transition, append the actual execution event. On completion, extract
allowlisted facts from the report using code; do not ask the LLM to invent a
summary status. Link failures and partial outputs too. The model may propose an
interpretation or next task separately. Validate that patch against the referenced
evidence before accepting it into current state, then atomically render PROJECT.md.

Never repeatedly summarize the previous summary as the only memory source. Rebuild
the concise view from canonical records and original evidence. Archive resolved or
superseded entries without deleting their lineage. A stale or missing artifact
invalidates dependent claims until checked; it does not erase the historical run.

## 7. Acceptance checks before calling this implemented

- An apparently successful but out-of-scope task is rejected or left as a proposal;
  successful execution alone cannot satisfy an unrelated mission requirement.
- Rephrasing the same calculation in a new chat reuses its compatible completed
  outputs or existing job; it does not dispatch an equivalent new calculation.
- Changing an export count reuses the ranking; changing a query invalidates only
  stages that depend on that query. Deliberate independent replicates remain possible.
- Failure of the intended workflow cannot silently change library coverage,
  protein assembly, requested outputs or the definition of completion.
- An explicit user revision takes effect and remains linked to the original
  mission; routine within-scope repairs do not trigger unnecessary approval loops.

- An old active constraint remains enforced after more than 100 later messages.
- The same task produces the same selected facts/constraints in a fresh chat and
  after restart, regardless of unrelated discussion added elsewhere in the project.
- Two constructs or sites under one target do not exchange references silently.
- Changing library/query versions prevents accidental reuse of incompatible ranks.
- Contradictory notes appear as unresolved; generated interpretations cannot promote
  themselves to verified observations or override sealed results.
- A queued/failed task cannot render as scientifically completed in PROJECT.md.
- Every reported metric is traceable to a task, artifact, field and unit where needed.
- Excess context produces explicit retrieval/omission metadata, never silent loss
  of mandatory constraints or truncation of literal molecular inputs.
- Restart, interrupted updates and concurrent chats do not lose revisions or submit
  duplicate tasks. A human edit does not silently override canonical state.
- "Continue this project" selects a ready task from the recorded plan or explains
  the concrete blocker; it does not construct a new target/story from chat fragments.

A document helps make work inspectable. Bounded retrieval, typed inputs, explicit
provenance and executable gates are what make its rules operational. This design
reduces opportunities for unsupported claims; it cannot guarantee a model will
never make an error. Live-provider behavior needs separate acceptance testing.
