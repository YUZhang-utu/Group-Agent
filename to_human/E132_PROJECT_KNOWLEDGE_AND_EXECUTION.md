# Target project knowledge and execution contract

Status: proposed architecture and operating specification, 2026-11-11.
This document does not mean the proposed dossier, context builder or execution
gates are already implemented. E131 remains the current runtime baseline.

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
