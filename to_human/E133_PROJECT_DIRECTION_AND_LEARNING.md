# Project direction and learning

E133 implements project continuity during ongoing work. Each target project keeps
its original mission, user-authorized revisions, requirements, feedback, evidence
and task ledger outside the conversation window. Completed work and corrections
inform later conversations; they are not limited to the initial project prompt.

## Use in the Agent

After pulling `feature/structure-guided-chat`, restart the Chat server and reload
the page. Open the intended target Project before creating a new conversation.
The sidebar's **Project direction & learning** opens the mission, deliverables,
constraints, proposed changes, learning records and task alignment/reuse records.
Existing projects receive revision 1 from their saved project metadata. Review
that initial objective if it was a generic workspace objective.

Example prompts, in the same target project:

```text
Record this project correction: keep the full supplied protein sequence and use
two identical protein copies. Propose the corresponding typed execution rule.
```

```text
Review what this project has already completed and its unresolved failures.
Use the saved structures and rankings to propose the next task toward our
candidate MOL2 deliverable. Do not repeat completed discovery or screening.
```

```text
The previous result used the wrong ligand reference. Record this correction
against task <task-id>, preserve the old result as historical evidence, and
explain which downstream tasks need new inputs before continuing.
```

```text
Import the complete ranking from task <task-id> into a candidate set named
screening-round-1. Export the first 100 molecules as MOL2 from this saved set.
```

Feedback quotes are taken from the actual current user message, with its message
and session identifiers. The model can propose structured mission changes; the
user applies the concrete proposal in the dialog. Routine within-scope work does
not require renewed mission approval. Conflicting feedback stays explicit rather
than silently replacing the original objective.

## What changes as work progresses

| Progress | Persisted record | Effect on later work |
| --- | --- | --- |
| User supplies a correction or preference | Exact sourced quote | Active corrections/preferences/constraints remain pinned beyond recent chat history |
| A task finishes | Actual execution state and allowlisted report fields | The next conversation can use saved inputs, results and artifact paths |
| A task fails or blocks | Specific error, task ID and report path | Diagnose or resume existing work rather than assume a scientific negative result |
| A dependency workflow advances | Completed nodes, pending dependencies and blockers | Context includes available next steps in the existing workflow |
| A result is reused | Source task and verified artifact inventory | Reused output is not counted as independent scientific evidence |
| User changes the mission | New revision preserving the original baseline | Future dispatch uses the new constraints; old tasks retain their execution revision |
| User accepts a requirement | Current revision, report evidence IDs and rationale | Task completion and user acceptance remain separate |

This is project-specific persistent memory and evidence reuse. No model weights
are trained. Improved provider reasoning or scientific accuracy has not been
measured. Feedback recording through natural language depends on the provider
calling the supplied tool; actual task transitions and selected report fields
are recorded automatically. User observations and model hypotheses are not
promoted to experimental facts.

## Durable document and bounded context

`<project_root>/PROJECT.md` is generated from canonical records in the Chat SQLite
database. It contains the original/current mission, requirements and acceptance
criteria, typed rules, feedback/evidence, task ledger and dependency workflows.
The database retains revisions and provenance; it is not a recursively rewritten
LLM summary. Preserve the Chat workspace database and project directories together
when backing up or moving projects.

Edit the mission through the dialog or an adopted proposal. If PROJECT.md is
manually changed, the generator reports `manual_edit_conflict` and preserves the
file. Canonical records still govern execution. Review and save those manual
edits separately before restoring a generated view; the agent does not silently
import arbitrary prose as authorized execution rules.

The model receives a project brief of at most 22,000 serialized JSON characters:
current mission, active user corrections/preferences/constraints, legacy pinned
constraint/decision notes, relevant/recent evidence, recent tasks and ready
workflow nodes. Full protein sequences are represented by length and hash in
this brief; exact inputs remain in canonical rules and scientific requests.
Optional evidence can be omitted with a count. Essential constraints that exceed
the budget cause an explicit error rather than silent loss. The broader dialogue
context has a 60,000-character budget, separate from system/tool schemas. These
are character limits, not token counts. Large old records remain available through
owned project search/report tools instead of being pasted into every turn.

Evidence records retain report hashes. Changed or unavailable reports are marked
stale. Report fields are recorded observations from software, not independent
validation of affinity, poses or biological activity.

## Execution checks and reuse boundaries

The Chat worker checks the concrete sealed plan before launching its subprocess.
Tasks link to current requirement IDs. Supported typed rules enforce permitted
or forbidden actions, literal AF3 sequences/ligands, identical protein copies,
model seeds, full-library versus budget-route restrictions and a pinned library
catalog version. A task started under an older mission revision is blocked on
resume until explicitly handled as a revised task; its sealed plan is not rewritten.

Free-text relevance, ligand-binding-site interpretation and arbitrary scientific
acceptance criteria are not mechanically proven by these checks. Explicit
requirement links plus typed rules cover the implemented constraints. Standalone
CLI execution retains its original checks; these new project checks apply to the
Chat worker. There is no claim of universal semantic drift detection.

Automatic content-based result reuse is limited to fingerprinted preparation and
search lanes: molecule preparation, literal AF3 input preparation, ligand-library
search and explicit target query/search requests. It includes inputs, parameters,
code, runtime/profile contents, package versions and library catalog/registry
versions. Completed artifact contents must still match their saved hashes.
Catalog versions identify the published library; this does not add a full shard
integrity scan on every chat turn. Keep library versions immutable.

Live database discovery, AF3 prediction execution and mutable candidate operations
retain their explicit result/reuse/resume mechanisms; E133 does not silently cache
those operations. A repeated eligible failed/in-progress task is blocked with its
existing ID. Reused output has a source-task link and launches no scientific
subprocess. Changing only the export count for an identical ligand-library search
blocks a new scan and directs the agent to import/export the saved ranking. This
guard does not automatically execute the replacement export or split every
combined workflow into cached stages.

If an independent repeat is intended, use this direct user command for the new
queued/blocked task, then explicitly resume it:

```text
/project_control {"operation":"replicate","task_id":"<new-task-id>","reason":"Independent repeat requested for comparison"}
```

The dialog lets the user archive resolved feedback and record requirement
acceptance with evidence. Acceptance applies only to that mission revision;
later revisions need review. The model cannot apply revisions or accept requirements
through its project_guidance tool.

## Validation scope

See [E133 protocol and results](../experiments/E133-project-guidance.md). Tests use
scripted providers, real local molecule preparation and a tiny real library for
search/ranking/MOL2 export. Live-provider behavior, GPU prediction, production
library throughput and browser visual acceptance are separate workstation checks.
