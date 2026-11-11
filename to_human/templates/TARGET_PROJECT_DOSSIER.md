# Target project dossier: <project name>

Template only. No fields below establish completed work. The proposed runtime
should render this view from versioned project records; that renderer is not yet
implemented. Use unknown/not_selected/not_run rather than inventing a value.

Primary purpose: keep every task aligned with the user's intended outcome.
Knowledge and progress records support this check; they must not gradually redefine
what the project was supposed to accomplish.

## Pinned mission and authorized changes

- Original user request / message ID: <verbatim request or durable source>.
- Original intended outcome: <what the user wants to achieve>.
- Current authorized mission revision: <revision and source user instruction>.
- In scope: <requested work and necessary prerequisites>.
- Out of scope: <explicit exclusions and unrelated activities>.
- Fixed requirements: <target, construct, assembly, chemistry, coverage, etc.>.
- Deliverables and definition of done: <observable outputs and checks>.
- Agent-proposed route: <method, separately editable within authorized scope>.

Preserve the original baseline. Only user-authorized changes replace goals, fixed
requirements or acceptance criteria; retain source messages and revision history.
Existing authorization counts. Routine method choices and repairs do not need
repeated approval. Agent suggestions remain proposals until authorized.

| Requirement ID | Required outcome / acceptance criterion | Current state | Evidence or blocker | User-authorized revision |
| --- | --- | --- | --- | --- |
| <R1> | <observable requirement> | outstanding | none | <source> |

## Direction check for the next task

- Requirement ID served: <ID>.
- Why this step advances that outcome: <concrete contribution or necessary dependency>.
- Expected deliverable or uncertainty resolved: <checkable result>.
- Scope and fixed requirements preserved: <checks>.
- Adjustment type: <within_scope_method_change, user_authorized_revision, or proposal_only>.
- Alignment verdict: <aligned, blocked, or out_of_scope>.
- Unrelated opportunities: <backlog references; do not silently execute>.

Repeat at dispatch, after results, on resumption and before selecting the next task.
Compare with the mission, not only the immediately preceding task. Failed execution
does not authorize weakening the goal or reporting a substitute deliverable as done.

## Identity and revision

| Field | Value |
| --- | --- |
| Project ID | <owned project ID> |
| Canonical state revision | <revision or manual draft> |
| Last updated | <timestamp> |
| Objective and requested deliverable | <user objective> |
| Target / accession / species | <values with provenance or unknown> |
| Construct / sequence artifact / hash | <exact artifact reference; do not paste long sequence here> |
| Assembly | <monomer, homodimer, heteromer or unresolved; explicit chain copies> |
| Active branch | <branch ID and purpose> |
| Current phase | not_started |

## Active constraints and decisions

Keep all applicable active constraints visible regardless of age. Link superseded
entries in history; do not simply delete them when a new message arrives.

| ID | Requirement / decision | Branch or input scope | Source message / evidence | Status / supersedes |
| --- | --- | --- | --- | --- |
| <ID> | <exact requirement> | <scope> | <source ID> | <active or proposed> |

## Evidence and reusable assets

Record only fields supported by the linked source. Input existence, payload
verification, runtime success and scientific validation are different checks.

| Evidence ID | Observation or asset | Task / artifact / field / hash | Verification level | Scope |
| --- | --- | --- | --- | --- |
| <ID> | <observation with units, or asset description> | <source> | <user_supplied, report_only, payload_checked, etc.> | <construct/site/library version> |

## Work plan and task ledger

This table is a current index. Full task cards and historical events remain in the
local records. Explicitly mark omitted stages not_applicable with a reason.

| Work item | Branch / dependencies | Input IDs and versions | Expected output / exit criterion | Execution state | Task / result reference |
| --- | --- | --- | --- | --- | --- |
| Establish target and construct | <branch> | <input IDs> | Identity and assembly recorded | not_started | none |
| Inspect reusable project work | <identity> | <prior task IDs> | Reuse and gaps recorded | not_started | none |
| Select reference/preparation/AF3 route | <inventory> | <selected assets> | Explicit route and choices | not_started | none |
| Execute selected scientific task | <ready dependencies> | <sealed inputs> | <task-specific completion/QC criteria> | not_started | none |
| Save/manage candidate collection, if applicable | <ranking> | <library/query/rank IDs> | Named immutable set and lineage | not_started | none |
| Deliver outputs and review next task | <checked outputs> | <result IDs> | Paths/counts/failures reported | not_started | none |

## Scientific interpretation and uncertainty

| Claim ID | Interpretation or hypothesis | Supporting evidence IDs | Limits / alternative explanation | Needed check |
| --- | --- | --- | --- | --- |
| <ID> | <explicitly tentative claim> | <IDs or none> | <limits> | <test or additional evidence> |

## Blockers and conflicting information

| Issue ID | Missing choice or conflict | Affected work items | Evidence on each side | Resolution required |
| --- | --- | --- | --- | --- |
| <ID> | <issue> | <IDs> | <sources> | <specific input/check> |

## Next ready action

- Work item and goal: <one concrete item, or blocked>.
- Mission requirement served and remaining gap: <requirement ID and contribution>.
- Reuse: <actual prior assets/tasks; no recomputation by default>.
- Required input and parameters: <typed values/references>.
- Completion evidence: <report fields/files/checks>.
- Stop condition: <specific missing information or failed check>.
- Scope of authorization: <requested execution, preparation only, or proposal>.

## Reuse and duplication check

Before creating a scientific task, compare against the actual project ledger.
Prompt wording and conversation IDs are not scientific work identities.

| Proposed work item | Equivalent prior tasks | Input/parameter/version compatibility | Artifact/checkpoint validity | Decision and reason |
| --- | --- | --- | --- | --- |
| <ID> | <IDs or no match> | <checked signature/differences> | <actual checks> | <reuse, wait, resume, recompute, or independent replicate> |

For changed inputs, list affected dependency stages and preserved reusable stages.
For a new calculation, explain why existing work cannot satisfy the requirement.
For uncertain submission, reconcile rather than submit again. A deliberate repeat
must retain the user's replicate intent.

## Task handoff card

Use for the selected work item; do not copy the entire project transcript.

```text
Project ID / revision:
Authorized mission revision / requirement IDs:
Contribution to the intended outcome:
Direction check / scope-change authorization if applicable:
Branch / work-item ID:
User request / source message:
Goal and deliverables:
Applicable constraint IDs:
Dependency task IDs and verified states:
Equivalent prior work / reuse decision and reason:
Changed inputs / affected stages, if any:
Input artifact IDs / hashes / versions:
Exact parameters and tool/runtime:
Output location:
Execution completion criteria:
Scientific quality criteria:
Unknowns / stop conditions:
Execution status / task ID:
Observed results with evidence pointers:
Proposed interpretation (separate from observations):
Next action or blocker:
```

## History and detailed evidence

- Canonical state / revision ledger: <local reference>.
- All task cards and execution reports: <local index>.
- Superseded decisions and resolved issues: <local history>.
- Searchable conversations: <project-scoped index>.
- Full structures, sequences, confidence arrays and candidate rankings:
  <artifact index; loaded only when relevant>.
