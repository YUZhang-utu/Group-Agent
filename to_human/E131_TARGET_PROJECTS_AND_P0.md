# Target projects, shared knowledge and durable scientific work

One Project represents a target. It contains multiple conversations and scientific
tasks, shared uploads, recorded decisions/constraints, candidate sets and durable
workflow goals. Switching a conversation does not discard the target's work.
Switching projects does not change the project of any queued or running task.

## Open the interface

Pull `feature/structure-guided-chat` and restart the existing Chat server with
the same storage directory, provider profiles and runtime configuration. Keep
`--allow-compute` enabled for searches/AF3 as before. Do not switch code versions
during a sealed running calculation. Project records reside in the local Chat
database; scientific inputs/outputs remain under the owned project directory.
Back up the entire storage directory, including its Chat/registry databases.

The left sidebar now has a target-project selector, New project, New conversation
and Project knowledge. Create a named project with target/accession and objective.
New conversations use the selected project. Tasks & results shows all tasks in
that project; individual chat transcripts remain separate. The knowledge panel
shows recorded notes, candidate collections and continuation status/errors.

Old sessions are migrated into their existing default project without moving
sealed inputs or outputs. This preserves provenance; it does not classify old
mixed-target chats automatically. Create separate projects for future target work.
Existing tasks are never rebound to another project.

## Prompt examples

Create/open the project first, then issue work in its conversation:

```text
Create a project named CIP2A with target Q8TCG1 and objective "Discover and
evaluate ligands for CIP2A". Keep all related work in this project.
```

```text
Show what this project has done, which tasks failed or remain blocked, and
which local structures and candidate sets can be reused. Search earlier project
conversations if the recent memory does not contain the answer. Do not rerun work.
Record this constraint: always preserve the complete supplied protein sequence.
```

For an attached single-molecule SMI/CSV/2D SDF (or a supplied SMILES), request a
durable goal instead of repeatedly prompting each asynchronous stage:

```text
In this project, create a durable research_workflow to prepare the attached
molecule with one conformer and seed 47, search every conformer in our configured
local library, retain the complete ranking as a named candidate set, and export
the first 1000 unique candidates as original and aligned MOL2. Continue successive
tasks automatically. If preparation yields multiple queries, show their actual
IDs and wait for my selection. Keep all results locally. Do not run docking.
```

```text
List this project's saved candidate sets. From set <SET_ID>, exclude molecules
labelled tested, retain at most one molecule per Murcko scaffold, and save the
result as "Untested scaffold representatives". Export its first 100 molecules
as MOL2 using the stored ranking and poses. Do not rescan the library.
```

Task IDs, set IDs, attachment IDs and query IDs must be returned by project tools,
not invented. The ordinary AF3 monomer/dimer/ligand prompts still run within the
selected project; this change does not introduce a new AF3 scientific model.
For a new experiment using a prior completed task's literal sequence/SMILES,
identify that task and ask to reuse its inputs. The workflow run decision records
the source task_id and carries its original user input into the new planning
request. Imported reports without recorded original input cannot use this path.

## Deterministic tool contracts

- `project`: create/open/list/memory/assets/note/search. Search covers persisted
  task requests, notes and conversation messages, with query/offset/limit.
- `molecule_prepare`: exactly one attachment_id or literal smiles. Optional
  conformers 1..100, seed, name and CSV column names. Default one conformer, seed47.
  At most 5000 records and 10000 requested conformers per task. Reports failures
  per record, duplicate chemistry, actual generated conformers and search queries.
- `candidate_sets`: import a completed ligand-library or budget ranking;
  list/show; union/intersection/difference; filter by rank limit or excluded labels;
  Gaussian score threshold; scaffold grouping; project labels; paged MOL2 export.
  Derived collections are immutable. Left/source order is retained; union appends
  unseen right members. There is no new scoring or cross-method score fusion.
- `research_workflow`: create/status/pause/resume/cancel/resolve/reconcile.
  Nodes specify id/tool/arguments/depends_on. The allowed queued tools are
  molecule_prepare, ligand_library, candidate_sets, workflow and target_search.
  `{{task:node}}`, `{{query:node}}`, `{{set:node}}` bind recorded task IDs, a sole
  prepared query or a produced candidate set. References must name dependencies.

The worker persists a dispatch marker before submitting a task. Following a crash
between submission and receipt, it blocks and lists candidate tasks for explicit
reconciliation; it never silently submits that step again. Pause/cancel affects
future dispatch, not a running scientific task. Resume a failed/interrupted task
explicitly using its existing ID, then resume its workflow. Resolve can replace
arguments of an unsubmitted node, for example to select one actual query ID.

## Chemistry and evidence limits

Original 3D coordinates are retained (one output conformer for each supplied 3D
record). Only generated conformers receive ETKDGv3 embedding and available MMFF
optimization. Optimization convergence is recorded. Salt fragments and unspecified
stereochemistry require explicit decisions; no automatic neutralization, tautomer
enumeration or silent component stripping occurs. A completed preparation task can
contain rejected records: inspect prepared_molecules, failed_molecules and
failures.csv. Zero or multiple usable queries block automatic query binding.

Library sets preserve molecule IDs, representative poses and source/library hashes.
MOL2 export verifies original records and chemical identity. Labels are recorded
statements, not measured activity. Scaffold grouping is exact nonstereochemical
Murcko grouping, not similarity clustering. Gaussian similarity is not affinity;
finite rigid-pose coverage does not establish active-compound recall. AF3 confidence
is not independent pose or biological validation.

The implementation is checked with local fixtures, real small-library scoring,
real worker subprocesses, HTTP and scripted-provider contracts. This does not
establish live-provider planning quality, GPU AF3 completion, production-library
speed/recall or completely unattended ambiguous target discovery. Windows browser
visual acceptance was blocked when Computer Use could not confidently establish
the current URL; no visual acceptance is claimed.
