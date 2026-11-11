# Target project dossier: <project name>

Template only. No fields below establish completed work. The proposed runtime
should render this view from versioned project records; that renderer is not yet
implemented. Use unknown/not_selected/not_run rather than inventing a value.

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
- Reuse: <actual prior assets/tasks; no recomputation by default>.
- Required input and parameters: <typed values/references>.
- Completion evidence: <report fields/files/checks>.
- Stop condition: <specific missing information or failed check>.
- Scope of authorization: <requested execution, preparation only, or proposal>.

## Task handoff card

Use for the selected work item; do not copy the entire project transcript.

```text
Project ID / revision:
Branch / work-item ID:
User request / source message:
Goal and deliverables:
Applicable constraint IDs:
Dependency task IDs and verified states:
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
