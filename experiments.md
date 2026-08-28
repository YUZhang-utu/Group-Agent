# Experiments

## E001 — Batch MOL2 registry bootstrap

Status: completed — confirmatory

### Hypothesis

A standard-library streaming parser can reliably register the initial library
without assuming one molecule per file, while terminal `_confN` naming is enough
to group the user's current conformers.

### Protocol

1. Create fixtures containing several MOL2 molecule blocks in one file.
2. Verify `name_conf0..N` records map to one molecule and distinct conformers.
3. Verify records without a suffix map to conformer zero.
4. Verify repeated imports are idempotent.
5. Verify two different structures claiming the same molecule/conformer index
   produce a visible conflict instead of silent overwrite.

### Acceptance criteria

- All unit tests pass.
- Import results are deterministic and queryable from SQLite.
- Existing source files are never modified.
- Invalid records are reported with file and record context.

### Classification

Confirmatory implementation validation.

### Results

- `python -m pytest -q`: 4 passed.
- Representative existing MOL2 dry-run: 1 file, 5 molecule records, 0 invalid.
- The source file was read-only and no production database was created by the
  dry-run.

The hypothesis was confirmed for the tested naming and MOL2 structures. A full
library import remains a separate operation because the user's actual library
location has not yet been selected.

### Production follow-up

The selected library at `6aa/` was subsequently validated and imported:

- Library ID: `LIB-AFA68EE6888C`
- Source files: 2
- Parsed/registered conformers: 299,999
- Unique molecules after `_confN` grouping: 100,000
- Invalid records/issues: 0/0
- Distribution: 99,999 molecules with three conformers; one molecule with
  conf1/conf2 only and no source conf0.

## E002 — Campaign and PDB target-structure selection

Status: completed — confirmatory implementation validation

### Hypothesis

An explicit, append-audited Campaign state machine can prevent target/structure
mix-ups while allowing RCSB candidates to be searched reproducibly and a human
to freeze one receptor structure with a scientific rationale before receptor
preparation begins.

### Protocol

1. Add versioned registry entities for Campaign, target, PDB candidates, and
   immutable decision events.
2. Enforce the transition order `draft -> target_set -> structures_review ->
   structure_selected`.
3. Query RCSB using a target identifier plus optional experimental-method and
   resolution constraints; preserve the query and returned metadata locally.
4. Require a non-empty rationale when selecting one candidate and prevent
   silent replacement of a frozen structure.
5. Expose the workflow through CLI commands and a machine-readable status
   summary.

### Acceptance criteria

- State-transition and invalid-transition tests pass.
- A mocked RCSB response can be registered deterministically without network
  access.
- Selection records PDB ID, rationale, timestamp, and decision event.
- Repeating candidate registration is idempotent.
- No receptor preparation or docking starts in E002.

### Classification

Confirmatory implementation validation. The eventual biological choice among
real PDB candidates remains a human-reviewed scientific decision.

### Results

- Added additive SQLite schema migration for Campaign, target, PDB candidate,
  and append-only decision-event records.
- Added CLI commands `create-campaign`, `set-target`, `search-pdb`,
  `select-pdb`, and `campaign-status`.
- Added an RCSB adapter with UniProt, method, and maximum-resolution constraints;
  it preserves resolution, method, deposition date, title, and ligand IDs.
- Seven full-project tests pass, including transition enforcement, idempotency,
  selection audit, and mocked RCSB metadata parsing.
- The production registry was migrated additively; Campaign tables remain empty
  pending a biological target definition.

The implementation hypothesis is confirmed. Scientific PDB selection is not
yet claimed: target identity and desired binding state must be defined first.

## E003 — Multi-user task-context foundation

Status: completed — confirmatory implementation validation

### Hypothesis

A mandatory `user -> project -> task -> campaign` context can isolate scientific
decisions and results between users while allowing models, docking engines, MD
engines, compound libraries, and compute caches to remain shared.

### Protocol

1. Add registry entities for users, projects, tasks, task access, and active
   task context.
2. Create an English-only task workspace with independent state, findings,
   experiment, benchmark, and research-log files.
3. Require every Campaign to belong to one task and validate user access before
   task activation or scientific work.
4. Add CLI commands to create, list, inspect, activate, and deactivate tasks.
5. Preserve the shared library registry and shared tool/model boundary.

### Acceptance criteria

- Two users cannot activate or inspect each other's private tasks.
- Tasks under different targets receive separate workspaces and continuity
  records.
- Campaign creation without a task is impossible through the public API/CLI.
- Existing molecule and conformer records remain unchanged after migration.
- The complete test suite passes.

### Classification

Confirmatory architecture and isolation validation.

### Results

- Added users, projects, tasks, task-access roles, and per-user active-task
  context to the shared registry.
- Added isolated English-only task workspaces with dedicated continuity,
  structure, docking, model, MD, data, and report paths.
- Campaign creation and every subsequent Campaign/PDB action now require the
  requesting user's active task to match the Campaign task.
- Added seven task-context CLI commands alongside the existing scientific CLI.
- Ten full-project tests pass, including cross-user access rejection, active
  task enforcement, separate workspaces, Campaign ownership, schema migration,
  and preservation of shared library records.
- The production registry migration preserved 100,000 molecules and 299,999
  conformers exactly.

The architecture hypothesis is confirmed for the local registry and CLI layer.
Authentication providers, collaborator-management UI, and platform-level role
administration remain later integration work.

## E004 — Project-centric optional workflows and provenance

Status: completed — confirmatory architecture validation

### Hypothesis

Replacing mandatory Task context with one active Project per user, while
recording only executed methods as generic Runs with immutable Artifact lineage,
will simplify group use without losing reproducibility or isolation.

### Protocol

1. Introduce a configured storage root with separate registry, shared, artifact,
   and per-user Project workspace boundaries.
2. Make Project the mandatory active scientific context; retain legacy Task
   records only for migration compatibility.
3. Make Campaign optional and bind it directly to a Project.
4. Register only methods that actually execute as workflow Runs. Record an
   explicit skip only when its rationale is scientifically meaningful.
5. Add Artifact, Run input/output, metric, and decision records with hashes and
   JSON manifests sufficient to reconstruct model, docking, rescoring, and MD
   provenance.
6. Reject paths that resolve outside the active Project or configured shared
   storage boundary.

### Acceptance criteria

- Users can activate only their own Projects and cannot read another user's
  private Project metadata through public APIs.
- Project workspaces are created under a fixed storage root with English-only
  continuity templates.
- A Project may contain any subset of method Runs; absent methods create no
  placeholder records.
- Run manifests register verified outputs and input/output lineage.
- Campaign/PDB actions require a matching active Project but no Task.
- Migration preserves the existing shared library, molecule, and conformer
  counts exactly.

### Classification

Confirmatory architecture and provenance validation.

### Results

- Added fixed storage-root initialization and isolated per-user Project
  workspaces with English-only continuity templates.
- Replaced Task with active Project in the public CLI and Campaign/PDB APIs.
  Legacy Task tables remain migration-only and are not exposed by the CLI.
- Added generic workflow Run, Artifact, input/output lineage, metrics, compute
  key, execution backend, and scientifically meaningful decision records.
- Added verified `result_manifest.json` import with SHA-256 calculation and
  strict Project path containment.
- Project summaries derive `executed_methods` only from actual Runs; an explicit
  skipped method does not appear as executed.
- Thirteen tests pass, including ownership isolation, optional-method behavior,
  manifest registration, metric persistence, path-escape rejection, Campaign
  Project binding, and legacy shared-library preservation.
- Production migration preserved exactly 100,000 molecules and 299,999
  conformers. No storage root, real user, or Project was invented.

The architecture hypothesis is confirmed. The next implementation boundary is
adapting individual model, docking, rescoring, and Slurm/MD executors to emit the
common result manifest.

## E005 — Cross-platform Slurm bridge

Status: completed — confirmatory architecture validation

### Hypothesis

A credential-free Slurm Bundle plus interchangeable OpenSSH and PuTTY command
transports can support Windows, macOS, and Linux clients while respecting daily
interactive cluster authentication and keeping the long-running Agent on the
desktop workstation.

### Protocol

1. Register immutable Cluster Profiles and per-user cluster usernames without
   passwords, MFA tokens, or private keys.
2. Generate deterministic, checksummed Run Bundles from structured Slurm
   resource specifications and registered remote roots.
3. Support OpenSSH on macOS/Linux/Windows and PuTTY Plink/PSCP on Windows.
4. Construct commands as argument arrays; never concatenate prompt text into a
   shell command.
5. Import a signed-by-content submission receipt containing Run, cluster,
   Slurm job, submitter, and Bundle identities.
6. Require interactive authentication when no reusable SSH session exists and
   never attempt to bypass daily cluster verification.

### Acceptance criteria

- Transport command tests cover OpenSSH and PuTTY without network access.
- Remote paths are restricted to registered POSIX roots.
- Generated `submit.slurm` contains only validated structured fields.
- Bundle manifests detect modification after export.
- Submission receipts cannot update another Project or Run.
- Existing Project, library, and provenance tests continue to pass.

### Classification

Confirmatory architecture and security validation.

### Results

- Added credential-free Cluster Profiles and per-user cluster username mappings.
- Added OpenSSH command generation for macOS, Linux, and Windows, plus PuTTY
  Plink/PSCP compatibility with saved-session connection sharing.
- Added structured Slurm resource validation, deterministic Bundle export,
  file hashes, path containment, and offline Bundle verification.
- Added Bundle registry and submission-receipt import; receipts must match the
  active Project, Run, Cluster, and exported Bundle hash.
- Added five public CLI commands for cluster registration, user linking, Bundle
  export/validation, and receipt import.
- Seventeen tests pass without network or cluster access. The production schema
  migration preserved 100,000 molecules and 299,999 conformers and created no
  fictitious Cluster Profiles, Bundles, or receipts.

The architecture hypothesis is confirmed at the control-plane boundary. Real
SSH authentication, transfer, `sbatch`, monitoring, and collection remain an
explicit integration test against the institution's cluster policy.

## E006 — Multi-structure pocket comparison and ligand similarity

Status: core implementation completed — chemistry deployment validation pending

### Hypothesis

A versioned comparison set containing multiple experimental or predicted
structures, combined with deterministic pocket superposition and persistent
ligand indices, can support fast interactive structure selection without
collapsing a target to one PDB ID prematurely.

### Protocol

1. Allow a Campaign to freeze multiple structure candidates into a named
   comparison set with one alignment reference.
2. Align structures on matched pocket C-alpha coordinates and report global
   pocket RMSD plus per-residue displacement and missing-residue coverage.
3. Generate a local PyMOL script that loads, aligns, groups, colors, and displays
   selected structures, pocket residues, and bound ligands.
4. Define a persistent RDKit 2D Morgan index and a two-stage 3D search using
   USRCAT prefiltering followed by exact shape reranking.
5. Keep RDKit, Gemmi, and PyMOL optional so the control plane remains usable on
   machines without scientific GUI/chemistry packages.

### Acceptance criteria

- A comparison set can contain multiple PDB/predicted candidates.
- Kabsch alignment recovers known rigid transforms and residue deviations.
- Missing pocket residues are explicitly reported rather than silently dropped.
- Generated PyMOL commands contain only validated IDs, paths, and residue lists.
- Similarity adapters fail with an actionable dependency message when RDKit is
  unavailable.
- Existing tests continue to pass.

### Classification

Confirmatory implementation and algorithm validation.

### Interim results

- Added comparison-set, comparison-member, and similarity-index registry
  entities without changing the existing single-candidate audit history.
- Added dependency-free Kabsch alignment, pocket RMSD, per-residue displacement,
  coverage, and missing-residue reporting.
- Added Project-scoped PyMOL script generation for multiple structures, chains,
  pockets, and organic ligands.
- Added lazy RDKit adapters for Morgan 2D search and USRCAT 3D descriptors plus
  explicit optional chemistry/structure dependency groups.
- Twenty-one tests pass, including rigid-transform recovery, local pocket
  deformation, missing residues, PyMOL path containment, and actionable missing
  dependency behavior.

Persistent high-throughput index construction and benchmarking remain pending
until canonical SMILES are attached to the library and the RDKit environment is
installed. RCSB mmCIF acquisition and ligand extraction also remain the next
implementation slice.

## E007 — Deployable target-structure and ligand-search toolchain

Status: core implementation completed — desktop integration validation pending

### Results

- Added Python 3.11 core and scientific-workstation Conda environments plus an
  environment diagnostic command.
- Added validated RCSB mmCIF acquisition, checksum metadata, AlphaFold 3 input
  generation, and external model-profile command construction.
- Added versioned disk-backed Morgan and USRCAT indices.
- Added source-control exclusions for all large/scientific/private artifacts.
- Twenty-seven dependency-independent tests pass. Real RDKit, PyMOL, and
  AlphaFold execution remains an explicit desktop/cluster integration test.

### Hypothesis

A dependency-layered desktop installation, project-scoped structure downloads,
validated AlphaFold job manifests, and disk-backed ligand indices can provide a
real testable workflow without committing scientific data, model weights, or
credentials to source control.

### Protocol

1. Define reproducible Conda environments for the core agent and chemistry /
   structure workstation, with an environment diagnostic command.
2. Download RCSB mmCIF files into the active Project, verify content, compute a
   SHA-256 checksum, and emit machine-readable acquisition metadata.
3. Generate AlphaFold 3 input JSON and AlphaFold 2 FASTA/Slurm artifacts from a
   validated, credential-free model profile; model weights remain external.
4. Build and query persistent Morgan 2D and USRCAT 3D indices with explicit
   format versions and deterministic ranking.
5. Exercise all dependency-independent paths in automated tests and document a
   desktop smoke test for dependency-backed paths.

### Acceptance criteria

- Invalid PDB identifiers, non-mmCIF responses, unsafe paths, and checksum
  mismatches fail before registration.
- AlphaFold configuration never stores weight files or credentials in Git.
- Generated jobs use only declared paths and structured resource settings.
- Similarity indices can be saved, loaded, and rejected on version mismatch.
- The repository excludes databases, molecular libraries, predictions, model
  weights, secrets, and user workspaces.
- Existing tests continue to pass.

### Classification

Confirmatory deployment and integration implementation.

## E008 — WEE1 real deployment validation

Status: in progress — confirmatory integration validation

### Hypothesis

The task-scoped agent deployed on university Linux storage can retrieve real
WEE1 structures, preserve Campaign provenance, and support human visual review
before a receptor is selected.

### Protocol

1. Create and activate one real User, Project, Task, and Campaign.
2. Bind reviewed human WEE1 metadata using UniProt accession P30291.
3. Search RCSB for X-ray candidates with resolution no worse than 3.0 angstrom.
4. Download selected candidates into the Project and open them from the
   terminal in PyMOL without changing Campaign selection state.
5. Define ATP-pocket residues, compare candidate pocket geometry and bound
   ligands, and record the final human selection rationale.

### Current results

- The real context was created and activated successfully.
- The WEE1 target was bound successfully.
- The RCSB search registered 25 candidates after one null-metadata parser fix.
- The full automated suite passes 29 tests.
- Visual PyMOL inspection and final structure selection remain pending.

### Classification

Confirmatory deployment validation. Structure choice remains a human-reviewed
scientific decision.
