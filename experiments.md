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
- Terminal-driven single-structure PyMOL inspection was validated successfully
  in the university environment.
- Multi-structure comparison, Codex natural-language orchestration, and final
  structure selection remain pending.

### Classification

Confirmatory deployment validation. Structure choice remains a human-reviewed
scientific decision.

## E009 — Campaign-aware multi-structure review CLI

Status: completed — confirmatory implementation validation

### Hypothesis

Exposing the existing comparison-set and PyMOL-review primitives through a
Campaign-aware CLI will make the real WEE1 structure-selection workflow
executable without copying internal database IDs or prematurely freezing one
receptor.

### Protocol

1. Let users name comparison members by PDB ID while resolving and validating
   their Campaign-scoped candidate IDs internally.
2. Require at least two PDB candidates, one member as the alignment reference,
   an explicit chain per member (default `A`), ATP-pocket residue numbers, and
   a scientific comparison rationale.
3. Add a read-only comparison-set status command with candidate/PDB metadata.
4. Generate a Project-scoped PyMOL script from downloaded
   `inputs/structures/<PDB>.cif` files and fail clearly when any file is absent.
5. Do not change Campaign state or call the irreversible final `select-pdb`
   checkpoint.

### Acceptance criteria

- PDB IDs outside the Campaign are rejected.
- The reference must be one of the selected members.
- Invalid chain specifications and missing structure files are rejected.
- The generated script loads every member, uses member-specific chains, shows
  the pocket and organic ligands, and aligns all mobile structures.
- Campaign remains in `structures_review` after comparison creation and review
  generation.
- The complete automated test suite passes.

### Classification

Confirmatory workflow-integration validation. The biological choice of WEE1
receptor remains human-reviewed.

### Results

- Added public PDB-ID resolution for Campaign-scoped comparison sets.
- Added `create-structure-comparison`, `structure-comparison-status`, and
  `prepare-pymol-review` CLI commands.
- Review generation validates downloaded mmCIF presence and writes only inside
  the supplied Project root.
- Comparison creation and review generation leave the Campaign in
  `structures_review`; only the separate `select-pdb` command can freeze it.
- The complete suite passes 31 tests, including foreign-PDB rejection,
  member-specific chains, missing-file failure, and Campaign-state preservation.

## E010 — Nonpolymer component classification

Status: completed — confirmatory implementation validation

### Hypothesis

Separating all PDB nonpolymer components from candidate binding ligands will
prevent common solvents, ions, buffers, and crystallization additives from
polluting Campaign ligand review while preserving the complete RCSB record.

### Protocol

1. Preserve every returned component ID in `nonpolymer_ids`.
2. Filter a versioned, conservative set of well-known solvent/additive/ion IDs
   from `ligand_ids`.
3. Record each exclusion and its reason in
   `excluded_nonpolymer_components` inside candidate metadata.
4. Keep unknown components as candidate ligands so an automated heuristic does
   not silently discard a genuine inhibitor or cofactor.
5. Add regression coverage for GOL, CL, NA, EDO, PO4, and MG, plus retention of
   a genuine ligand ID.

### Acceptance criteria

- Known non-ligand components remain visible in full metadata but not in
  `ligand_ids`.
- Unknown and drug-like component IDs remain in `ligand_ids`.
- Null RCSB nonpolymer lists remain supported.
- Existing Campaign storage requires no destructive schema migration.

### Results

- RCSB parsing now emits `nonpolymer_ids`, filtered `ligand_ids`, and
  `excluded_nonpolymer_components` with explicit reasons.
- GOL, CL, NA, EDO, PO4, and MG are retained as nonpolymers but excluded from
  candidate binding ligands; ATP and unknown component IDs remain candidates.
- Repeating `search-pdb` now refreshes metadata for existing Campaign
  candidates while its `inserted` count continues to report only new rows.
- `campaign-status` exposes both the filtered and complete component views.
- Full suite: 32 tests passed.

## E011 — AI recommendation and prompt orchestration boundary

Status: completed — confirmatory architecture validation

### Hypothesis

A provider-neutral, schema-validated recommendation layer can make prompts the
primary Project interface while preserving deterministic execution, scientific
provenance, user approval checkpoints, and the exclusion of experimental raw
data from model context.

### Protocol

1. Register immutable AI review requests owned by one active Project and
   optional Campaign.
2. Store task type, subject, prompt version, rendered prompt, evidence packet,
   declared data classes, and model-facing privacy policy.
3. Reject requests declaring experimental raw data or experimental
   measurements as model inputs.
4. Import a structured recommendation containing action, rationale, confidence,
   evidence references, uncertainties, and whether human review is required.
5. Record provider/model identity and preserve the model response without
   directly mutating Campaign or Run state.
6. Expose create/show/import operations through CLI commands so Codex or any
   later provider adapter can use the same audited boundary.

### Acceptance criteria

- Cross-Project access is rejected.
- Forbidden experimental data classes are rejected before request creation.
- Invalid confidence, missing rationale, and unreferenced evidence are rejected.
- Importing a recommendation does not select a PDB or execute a Run.
- Requests and recommendations round-trip through JSON and the complete suite
  remains green.

### Classification

Confirmatory architecture and privacy-boundary validation.

### Results

- Added immutable `ai_review_request` and `ai_recommendation` registry entities.
- Added provider-neutral create/import/status APIs and CLI commands.
- Requests preserve rendered prompts, versioned evidence, declared data
  classes, and a machine-readable privacy policy.
- Raw experimental data and experimental measurements are rejected as model
  context; unknown data classes are also rejected.
- Recommendations require valid confidence, rationale, evidence references,
  uncertainty, provider/model identity, and a human-review flag.
- Importing a recommendation leaves Campaign and Run state unchanged.
- Full suite: 35 tests passed.

## E012 — Predicted target-structure acquisition and backend deployment

Status: completed — confirmatory deployment-adapter validation

### Hypothesis

A provider-neutral structure-prediction adapter covering AlphaFold DB,
AlphaFold 3, Boltz-2, and Chai-1 can produce versioned WEE1 structure candidates
without coupling the target-preparation workflow to one transient
state-of-the-art model.

### Protocol

1. Add validated AlphaFold DB mmCIF acquisition by UniProt accession with hash
   and source metadata.
2. Generalize local model profiles to AlphaFold 2/3, Boltz-2, and Chai-1 while
   retaining explicit license acknowledgement where required.
3. Generate backend-native, Project-scoped inputs from the same reviewed amino
   acid sequence.
4. Construct commands as argument arrays, keeping model caches, weights,
   databases, and credentials outside Git and the LLM context.
5. Record backend/model identity so outputs can later enter the common target
   structure registry as immutable candidates.

### Acceptance criteria

- Invalid UniProt IDs and invalid model profiles fail before execution.
- AFDB responses are verified as the requested mmCIF before writing.
- AF3 JSON, Boltz YAML, and Chai FASTA inputs are deterministic and
  Project-scoped.
- Commands match each backend's official CLI boundary and never embed secrets.
- Existing AlphaFold 2/3 profiles remain compatible.
- Full tests pass without requiring GPU models or network access.

### Classification

Confirmatory deployment-adapter validation; real WEE1 model generation remains
a workstation/cluster integration run.

### Results

- Added verified AlphaFold DB mmCIF acquisition by UniProt accession.
- Added a unified Project-scoped input layer for AlphaFold 2/3, Boltz-2, and
  Chai-1.
- Added profile validation and argument-array command generation for Boltz-2
  and Chai-1 while preserving existing AlphaFold profiles.
- Added example Boltz-2 and Chai-1 deployment profiles with external model/cache
  paths and no embedded credentials.
- Full offline suite: 38 tests passed.
- Real WEE1 inference remains pending on the university GPU environment and is
  not claimed by this implementation result.

## E013 — Predicted-structure Campaign registration and receptor lock

Status: implemented — workstation import pending

### Hypothesis

An immutable, model-agnostic prediction-result manifest can register a real
AlphaFold 3 output beside experimental PDB candidates without treating either
source as automatically authoritative, while a single receptor-selection lock
preserves the exact file hash, construct, model provenance, and human rationale.

### Protocol

1. Inspect a completed prediction output directory and identify its final mmCIF,
   ranking table, confidence summary, input, and content hashes.
2. Register predicted structures in a Campaign-owned table separate from RCSB
   candidates so PDB semantics are not overloaded with synthetic identifiers.
3. Expose experimental and predicted candidates together in Campaign status.
4. Add a generic receptor-selection operation accepting either candidate type
   and write an immutable version-lock snapshot before changing Campaign state.
5. Preserve the existing `select-pdb` interface as a compatibility wrapper.
6. Keep WEE1 as an integration example only; no target name, PDB, path, or model
   backend may be hard-coded into the registry API.

### Acceptance criteria

- Modified or missing structure files are rejected during registration.
- Duplicate prediction imports are idempotent by Campaign and structure hash.
- A selected predicted structure records backend/version, construct, score,
  source paths, and checksums in the immutable lock.
- Campaign selection still requires an explicit scientific rationale and human
  action; importing an AF3 result cannot select it automatically.
- Existing PDB Campaign behavior and the complete offline suite remain green.

### Classification

Confirmatory integration. The real WEE1 AF3 run and PyMOL RMSD observations are
deployment evidence supplied after E012, not retrospective E012 test results.

### Results

- Added hashed AlphaFold 3 output inspection for the final mmCIF, ranking CSV,
  confidence summary, input JSON, and runtime metadata.
- Added Campaign-owned predicted candidates without overloading PDB identifiers.
- Added generic experimental/predicted receptor selection with an immutable
  snapshot lock and preserved `select-pdb` compatibility.
- Added `import-alphafold3-result` and `select-receptor` CLI operations.
- Offline suite: 41 tests passed.
- Real WEE1 deployment evidence: end-to-end AF3 inference succeeded on an RTX
  5090; the best experimental/PDB comparison gave approximately 1.10 Å C-alpha
  RMSD without outlier rejection and 0.37 Å after rejection. Import into the
  production Campaign remains pending because the production registry and AF3
  files reside on the university Linux workstation.

## E014 — Unified receptor ensemble and ligand-pocket review

Status: implemented — workstation validation pending

### Hypothesis

Using one experimental structure only as an alignment reference while retaining
all experimental and predicted candidates will expose receptor-state and ligand-
pose diversity before selection, without implying that the reference is best.

### Protocol

1. Create a Campaign ensemble containing all PDB and predicted candidates, with
   one explicit reference candidate and per-candidate chain assignments.
2. Preserve filtered ligand IDs from RCSB; predicted apo candidates have no
   invented ligand.
3. Generate a PyMOL review that aligns every protein to the reference, displays
   all retained ligands separately, and selects each ligand's local protein
   neighborhood.
4. Keep canonical pocket residues as an optional shared comparison selection.
5. Leave the Campaign in `structures_review`; ensemble creation and visualization
   must never select or lock a receptor.

### Acceptance criteria

- `8BJU` can be the reference while every Campaign PDB and AF3 candidate remains
  visible as a peer.
- Foreign candidate IDs, unsafe chains, and missing files are rejected.
- Excluded solvent/ion component IDs do not reappear as ligands.
- Generated scripts distinguish experimental ligands from predicted apo models.
- Offline tests pass and no target-specific identifier is hard-coded.

### Results

- Added all-Campaign receptor ensembles containing every experimental PDB and
  every registered prediction, with one explicit alignment reference.
- Added reference-neutral status and Project-scoped PyMOL generation.
- Generated scripts use `cealign`, display only classified ligand IDs, build
  per-structure 6 Å ligand pockets, and retain optional canonical pocket
  selections.
- Predicted construct names ending in `_START_END` automatically map canonical
  residue numbering to model-local numbering; WEE1 299–569 therefore uses an
  offset of 298.
- Ensemble creation leaves Campaign state at `structures_review`.
- Added rationale-required ensemble exclusions so UniProt-matched fragments or
  off-domain complexes remain in Campaign provenance without entering an
  incompatible receptor alignment. The motivating case is 9TG7: WEE1 is a
  12-residue degron peptide rather than the kinase domain.
- Offline suite: 42 tests passed. Real 8BJU/all-WEE1 visualization remains the
  workstation validation step.

## E015 — LLM receptor-eligibility recommendation interface

Status: implemented — workstation validation pending

### Hypothesis

A target-agnostic evidence packet and strict per-candidate response schema can
let an LLM identify short fragments, off-domain complexes, or wrong chains while
keeping final exclusion under deterministic validation and human control.

### Protocol

1. Build an AI request from Campaign target metadata, receptor-purpose/domain
   requirements, public PDB metadata, prediction provenance, and optional
   computed coverage summaries.
2. Require one `include`, `exclude`, or `manual_review` decision per candidate,
   each with rationale and evidence references.
3. Reject unknown/missing candidates, invalid verdicts, missing evidence, and
   any response that disables human review.
4. Store provider/model/version and the complete response using the existing
   immutable AI recommendation registry.
5. Do not mutate an ensemble or Campaign; accepted exclusions remain a later
   explicit human action.

### Acceptance criteria

- The interface contains no WEE1- or 9TG7-specific rule.
- Experimental raw data and measurements remain forbidden from model context.
- A 9TG7-like evidence record can be recommended for exclusion with cited
  chain-length/domain evidence.
- Partial or hallucinated candidate lists are rejected.
- Existing generic AI review behavior remains compatible.

### Results

- Added a provider-neutral request builder for Campaign target requirements,
  public PDB metadata, prediction provenance, and optional computed evidence.
- Every candidate must occur exactly once in the response. Verdicts are limited
  to `include`, `exclude`, or `manual_review`, evidence citations are checked,
  and human review cannot be disabled.
- The dedicated CLI emits the complete prompt/evidence packet for a later LLM
  provider adapter. Recommendations never mutate Campaign state or ensembles.
- Confirmatory 8BJU/9TG7-like tests use domain-purpose evidence without
  hard-coding either identifier. Full offline suite: 44 tests passed.

## E016 — Campaign ligand registry, comparison, and query lock

Status: implemented — workstation chemistry validation pending

### Hypothesis

Representing each retained co-crystal ligand as a provenance-rich, standardized
Campaign entity will prevent CCD identity, crystal instance, and normalized
chemical identity from being conflated, while enabling an auditable human query
selection for downstream similarity search.

### Protocol

1. Register ligand instances from an explicit manifest containing PDB/CCD,
   chain/residue identity, source structure, standardized SMILES, and optional
   crystal conformer artifact/checksum.
2. Reject components not present in the parent candidate's filtered ligand list;
   never revive excluded solvents, salts, buffers, or ions.
3. Compare registered ligands using Morgan/Tanimoto and, where 3D conformers are
   available, USRCAT; represent missing 3D data explicitly rather than inventing
   conformers.
4. Lock one query ligand with a required rationale and immutable chemical/
   provenance snapshot. Registration and comparison must not select it.
5. Expose provider-neutral AI evidence for later recommendation without allowing
   an LLM to create the lock.

### Acceptance criteria

- Multiple instances of the same CCD ligand remain distinguishable.
- Filtered-out PDB components cannot be registered as Campaign ligands.
- 2D and 3D comparisons report their methods and parameters independently.
- Query selection is user-authorized, rationale-required, and immutable.
- Existing receptor and AI workflows remain green.

### Results

- Added provenance-rich Campaign ligand instances with PDB/CCD, chain, residue,
  altloc, standardized SMILES, optional source checksum, and optional USRCAT.
- Registration is constrained by each PDB candidate's retained ligand list, so
  excluded crystallization components cannot re-enter through this interface.
- Added pairwise Morgan and available-USRCAT comparison, an immutable query
  ligand snapshot lock, and a provider-neutral LLM query recommendation request.
- LLM assessments must cover every ligand and cannot create the lock.
- Offline core suite: 49 tests passed; real RDKit comparison remains a
  workstation chemistry-environment validation.

## E017 — Hierarchical 2D/3D macrocycle-library search

Status: implemented core — production index validation pending

### Hypothesis

A persistent Morgan/Tanimoto first stage followed by conformer-aware USRCAT
reranking can retrieve chemically and shape-relevant macrocycles reproducibly
without an all-against-all exact 3D alignment.

### Protocol

1. Bind ready Morgan and USRCAT indices to an explicitly locked query ligand and
   registered library.
2. Retrieve a configurable 2D candidate pool, then rerank available conformers
   by 3D similarity while retaining molecule/conformer identity.
3. Record index versions, parameters, cutoffs, limits, query snapshot, and ranked
   results in an immutable search record.
4. Keep 2D-only hits when 3D data are unavailable, labeling the missing stage.
5. Produce deterministic JSON suitable for later LLM explanation and human
   review, never automatic docking submission.

### Acceptance criteria

- Repeated searches with identical query/index/parameters return identical order.
- Molecule-level results retain the best conformer and both 2D/3D scores.
- Search rejects cross-library indices and unlocked/unowned query ligands.
- No experimental measurements enter LLM evidence.
- Offline tests cover query locking, filtering, ranking, and provenance.

### Results

- Added deterministic Morgan-pool/USRCAT reranking with best-conformer retention
  and explicit 2D-only fallback.
- Added immutable Campaign search records containing query snapshot, library,
  index paths and hashes, parameters, stages completed, and ranked results.
- The production library still needs chemistry-side standardized SMILES and
  conformer descriptors to build its real indices; existing MOL2 registration
  alone does not fabricate those representations.

## E018 — Authoritative ligand preparation and production index build

Status: protocol locked — implementation and workstation validation pending

### Hypothesis

Campaign mmCIF coordinates can be joined to cached RCSB CCD definitions by
atom name, preserving authoritative CCD bond orders, and the registered MOL2
records can be streamed through one versioned RDKit standardization policy to
produce reproducible molecule-level Morgan and conformer-level USRCAT indices
without inventing chemistry or coordinates.

### Protocol

1. Enumerate only retained Campaign ligand instances from `_atom_site`, keeping
   model, author chain/residue/insertion code, and explicit alternate location.
2. Fetch each CCD definition from the official RCSB ligand endpoint into a
   content-addressed cache; reject component-ID mismatches and missing bonds.
3. Join crystal atom names to CCD atoms, use only CCD bond orders, standardize
   the parent with an explicit RDKit policy/version, and write one SDF plus a
   checksummed manifest record per valid crystal conformer.
4. Compute USRCAT only when the joined molecule has a valid 3D conformer;
   report skipped/invalid instances instead of generating coordinates.
5. Stream every registered MOL2 conformer through RDKit, require conformers of
   one registered molecule to agree on standardized parent SMILES, and build a
   molecule Morgan index, conformer USRCAT index, conformer-to-molecule map,
   and checksummed manifest.
6. Make preparation and index construction resumable CLI operations whose
   outputs remain outside Git and whose manifests capture input/output hashes,
   software versions, parameters, counts, and failures.

### Acceptance criteria

- Retained-component filtering and full crystal instance identity are tested.
- CCD bond orders, never distance guesses, define ligand connectivity.
- Missing CCD atoms/bonds, invalid 3D records, and inconsistent MOL2 parents
  fail or skip visibly according to the locked policy.
- Repeated builds have deterministic ordering and manifests.
- The existing 49-test suite remains green without requiring network, RDKit,
  Gemmi, or production data; dependency-required paths fail clearly.

### Classification

Confirmatory implementation validation. The real WEE1 extraction and full
`LIB-AFA68EE6888C` build remain production workstation runs.

## E019 — Million-scale FAISS/Gaussian retrieval prototype

Status: protocol locked — implementation pending

### Hypothesis

An in-memory FAISS IVF-PQ index over USRCAT can provide a broad, alignment-free
L1 recall set from a billion conformers, after which local-NVMe shape metadata
and joint Gaussian shape/color optimization can recover accurate rankings.
Protein-interaction anchors, projected polar points, and a pocket exclusion
grid advance only if controlled ablations justify their complexity.

### Protocol

1. Calibrate the immutable library conformer protocol on Platinum: symmetry-
   corrected minimum heavy-atom RMSD to the bound pose, overall and stratified
   by rotatable bonds, heavy atoms, ring systems, molecular weight, and a
   query-chemotype-matched subset. Also report library conformer-count
   distributions within the same strata.
2. Audit every crystal query before feature extraction: resolution, occupancy,
   atom B factors, available RSCC/EDIA, CCD bond orders, added hydrogens, and
   protonation compatibility with the library. Unsupported polar atoms cannot
   become anchors.
3. In one sequential, rate-limited pass over read-only NFS shards, generate all
   local artifacts together: int16 heavy-atom coordinates at 0.01 angstrom,
   int16 feature coordinates and types, USRCAT float32 vectors, binary metadata,
   heavy-atom bounding boxes, and PMI axes. Each shard has an atomic done marker
   and checksums. Never fetch coordinates from NFS during online search.
   Store heavy-atom count and six parent-feature counts (HBD/HBA/cation/anion/
   hydrophobe/aromatic as uint8) in the L1.5 record; alternative projected
   directions do not count as separate functional features.
   Artifacts are immutable append-only shards registered by a small catalog;
   adding molecules creates new shards and stable global int64 conformer IDs
   without rewriting old data. FAISS uses `add_with_ids` for new shards. Track
   transformed-vector centroid/variance and coarse-list occupancy drift; retrain
   a new versioned IVF-PQ index only when measured recall or drift crosses a
   preregistered threshold, while retaining the previous index for rollback.
4. Preserve every supplied library conformer; do not RMSD-deduplicate or cap
   conformer counts. Prefer SDF where equivalent SDF and MOL2 sources exist.
   Parse independent shards concurrently and throttle aggregate NFS reads to an
   operator-configured fraction of measured bandwidth.
5. Separate graph chemistry from coordinate ingestion. Perform the minimum
   documented RDKit sanitization needed for aromaticity, formal charge,
   hybridization, hydrogen addition, and pharmacophore typing once per unique
   molecular graph; reuse it across conformers. A `sanitize=False` coordinate
   fast path cannot itself define chemical features.
   First audit whether conformers are contiguous and which source ID reliably
   defines molecule identity. The reusable template contains graph atom roles
   and projection rules, while projection coordinates are recomputed for every
   conformer from that conformer's geometry.
6. Use six feature types: HBD, HBA, cation, anion, hydrophobe, and aromatic-ring
   center. Use two query weight classes only: anchor=2.0, ordinary=1.0;
   solvent-exposed polar features are removed. Place polar anchor color points
   at validated ideal/protein-side hydrogen-bond projection points.
7. Define color Gaussian overlap with sigma=1.0 angstrom, same-type matches
   only, and skip pairs beyond 4.5 angstrom. Report both standard ColorTanimoto
   and query-biased ColorTversky (alpha=0.95, beta=0.05). Do not reuse standard
   ROCS Combo thresholds for weighted scores; calibrate cutoffs empirically.
8. After all USRCAT vectors are durable, train FAISS IVF-PQ on a fixed random 1%
   sample and add all vectors in resumable batches. Store per-dimension mean/std
   from that sample and z-score both database and query vectors identically;
   apply optional HBD/HBA block weights after z-scoring so normalization does
   not cancel them. Compare d60/m20, d60/m30, and zero-padded d64/m32 at nbits=8
   and initial nlist=32768. Retain the 240-GB raw float32 vectors for retraining.
   Explicitly load without mmap, verify direct_map is disabled, and warm up.
   After add/write, restart a clean process before steady-state memory/query
   benchmarks so build buffers and allocator fragmentation are released.
9. L1 retrieves a configurable 300,000--1,000,000 conformers. L1.5 uses stored
   bounding boxes, volume/extent ratios, and PMI eigenvalue ratios to reduce
   the set to an initial target near 100,000 without reading NFS. Under an
   injective parent-feature coverage rule, first reject candidates whose six-
   type feature counts cannot cover the required query anchors.
10. Generate rigid seeds without a triplet index: four proper PMI-axis sign
   assignments plus anchor-centroid/axis alignments with six sampled rotations.
   Deduplicate transforms and jointly optimize rotation/translation for shape
   and color from each seed. Prototype 10 seeds and 12--15 L-BFGS iterations,
   but treat the claimed 2--5 minute latency as a benchmark target, not fact.
11. Apply a 0.5-angstrom protein excluded-volume grid after overlay. For the
   final 10,000 candidates only, optionally relax one or two terminal rotatable
   bonds connected to polar groups and recompute the joint score.
12. Prototype on a one-million-molecule subset with Python orchestration and
   vectorized/compiled primitives. Measure before fixing the C++/Rust/CUDA
   production representation.
13. Generate directional interaction-site projections deterministically. HBD
   sites extend 1.9 angstrom from each explicit H along D--H. HBA sites extend
   2.9 angstrom from the acceptor along enumerated lone-pair directions: two
   in-plane +/-60 degree sites for carbonyl-like sp2 oxygen, two tetrahedral
   sites for sp3 oxygen, and one in-plane site for pyridine/imidazole-like sp2
   nitrogen. Ions remain at atom/group centroids. Aromatics initially use ring
   centers; +/-3.5 angstrom normal projections are a separate ablation.
14. Use the same versioned protonation, tautomer, hydrogen-addition, and local
    geometry policy for CCD-derived queries and library molecules. Record it in
    every feature/index manifest. Enumerate ambiguous donor-H/lone-pair
    directions explicitly rather than choosing an undocumented orientation.
15. Treat multiple projection points as one functional feature in weighting.
    Split its functional weight equally across its `n` alternative projected
    directions (`W/n` per point), and record the parent feature ID so reporting
    and ablations remain group-aware. Because this mixture convention changes
    self-overlap relative to a single point, calibrate all Tanimoto/Tversky
    thresholds empirically and never compare them with ROCS values.
16. Store ShapeTanimoto, query-biased ShapeTversky, ColorTanimoto, and
    `color_tversky_qbiased` separately. Define the Tversky operand direction in
    the schema. Call the experimental sum `anchor_combo`, not ROCS
    ComboTanimoto, and mark it incomparable to ROCS thresholds.
17. Calibrate retain/priority/high-similarity cutoffs on 3--5 development
    targets with multiple known actives. Use one cocrystal query per target,
    held-out actives, and at least 100,000 fixed random molecules. Initialize
    retention at the random-score 99.9th percentile and report fixed-N hits;
    freeze thresholds before evaluating held-out targets. Bind every cutoff to
    the exact feature/projection/weight/protonation/score configuration.
18. Write one immutable result directory per query/configuration hash. Store
    molecule ID, conformer ID, component scores, anchor RMSD, seed provenance,
    and final 4x4 rigid transform for every ranked result. Export overlaid SDF
    only for the top 1,000 by default.

### Required ablations and acceptance criteria

- FAISS candidate recall@N against exhaustive joint scoring on a tractable
  subset; report recall@100, @1000, and @10000 across nlist/nprobe/PQ settings.
- Artifact bytes/conformer, NFS ingest throughput, parsing throughput, FAISS
  train/add time, cold/warm latency, L1/L1.5 candidate counts, mmap page faults,
  and every L2/L3/L4 stage latency percentile.
- If cold-cache local record reads dominate, sort candidates by file offset for
  batched near-sequential reads and restore result order by conformer ID.
- Shape-only, shape+unweighted-color, shape+anchors, projected-vs-atom-centered
  anchors, Tanimoto-vs-query-Tversky, and terminal-torsion-refinement ablations.
- DUD-E and LIT-PCBA retrospective evaluation with EF1%, BEDROC, hit rate at
  fixed N (especially 1,000 and 10,000), per-target results, and uncertainty.
  Interpret LIT-PCBA as relative method-vs-baseline evidence, not absolute
  binding truth. DUD-E is secondary because of known analogue/decoy biases.
- Mandatory ordered baselines: Morgan 2D; USRCAT; aligned shape-only Gaussian;
  shape plus unweighted atom-centered color; and the full anchor-weighted,
  projected-color scheme. The 3->4 and 4->5 deltas determine whether color and
  then anchors/projections justify their complexity.
- The complex method advances to a native billion-scale engine only if index
  recall is acceptable and color/anchor variants materially improve early
  enrichment over shape-only and USRCAT baselines.
- Molecular overlap uses one documented approximation order consistently for
  cross-overlap and both self-overlaps.

### Classification

Confirmatory benchmark and architecture-selection experiment. All parameter
tuning discovered after protocol lock is exploratory until repeated on a held-
out benchmark split.

### Exploratory laptop smoke results (2026-09-02)

### Protocol amendment: anchors are reranking evidence only (2026-09-02)

The earlier hard anchor-count L1.5 gate and zero-weight solvent-exposed feature
rule are withdrawn. Observed protein--ligand hydrogen bonds are structural
facts, but their biological necessity is unknown; anchor-derived information
must not cause irreversible candidate deletion.

1. L1 remains anchor-independent USRCAT/FAISS with a deliberately broad pool.
   L1.5 may hard-filter only anchor-independent geometry such as very loose
   volume and PMI aspect-ratio bounds. Feature counts are annotations or
   ranking signals, not admission criteria.
2. Default query weights become anchor=1.5, ordinary=1.0, and
   solvent-exposed=0.5. No chemically valid polar feature receives zero weight.
3. L2 stores the transform, shape scores, unweighted atom-centered color
   scores, anchored color, atom-centered and projected color, and one separate
   overlap value for every query anchor for every retained conformer.
4. Atom-centered and projected representations remain independently normalized
   scoring views. Projected color cannot replace atom-centered color or gate a
   candidate. Initial ranking uses atom-centered color; max-color and
   projected-only remain reversible evaluation views.
5. L2 output count is configured from downstream capacity, initially 2--5x the
   next stage's consumption. No universal fixed cutoff is chemical truth.
6. Compare anchored and unweighted rankings with Spearman correlation and
   top-1,000 overlap. Above 90% overlap motivates simplification; below 60%
   triggers discordant-hit review before anchor ranking is promoted.

This amendment supersedes E019 steps 5 and 9 only where they specify zero
weights or anchor-count filtering. The original text remains for auditability.

### Result-schema amendment: anchor identity and objective-specific poses

1. Every query output directory contains `query_manifest.json`. Each anchor has
   a stable ID, ligand atom indices, feature type, atom center, all projected
   point coordinates, evidence fields (interaction type, distance, angle, and
   burial), and the weight used for the materialized anchored-score column.
2. A candidate contains one named pose per optimized objective, for example
   `shape_only`, `atomcentered_joint`, and `projected_joint`. Every pose owns a
   separate 4x4 transform and a full set of shape/color scores evaluated under
   exactly that transform. Scores from different optimal transforms cannot be
   silently combined.
3. `per_anchor_overlap_raw[i]` stores the unnormalized query-anchor/candidate
   cross-overlap and `per_anchor_self_overlap_raw[i]` stores the corresponding
   query self-overlap. Position `i` maps to the ordered anchor definitions in
   `query_manifest.json`; reporting and reranking expose stable anchor IDs.
4. The 1.5/1.0/0.5 weights are a recorded materialization preset, not a search
   invariant. Raw component data support later weight/normalization changes
   without rigid realignment. Any derived anchored column records its weights.

The first resource-scaled test used synthetic, correlated, scale-heterogeneous
60-D vectors. It validates the benchmark machinery and memory safety only; it
does not establish performance on real USRCAT or at billion scale.

- Laptop: 32 logical CPUs, 31.05 GiB RAM, 16.97 GiB available at start.
- 200,000 vectors, 20 queries, IVF1024-equivalent small setting (nlist=512,
  nprobe=32), exact top-100 covered by returned top-1000: all three PQ variants
  reached 100% candidate recall, so the setting was not discriminative.
- Strict run: 1,000,000 vectors, 50 queries, nlist=1024, nprobe=8, returned
  top-100 evaluated against exact top-100.
- d60/m20: recall 0.9010, 28.32 serialized bytes/vector.
- d60/m30: recall 0.9464, 38.32 serialized bytes/vector.
- padded-d64/m32: recall 0.9474, 40.34 serialized bytes/vector.
- The padded variant gained only 0.001 absolute recall over m30 while adding
  about 2.02 bytes/vector in this synthetic run. m30 is the provisional default,
  but no production choice is made until real-vector and nprobe/top-N sweeps.
- Strict-run wall time was 16.45 seconds; peak observed process RSS among the
  reported configurations was 0.864 GiB. These timings include only synthetic
  vector generation, exact neighbor truth, and FAISS—not chemistry or L2.

### Real 100K-molecule/299,999-conformer results (2026-09-02)

- Query: WEE1 8BJU QT9, model 1, author chain A, residue 601. The 1.53-angstrom
  structure's validation report gives ligand RSCC 0.933, EDIA mean 0.951,
  occupancy 1.000, 42/42 density-supported atoms, and mean B factor 26.74
  square angstrom. CCD supplied bond orders; extraction produced 42 heavy atoms.
- Full library RDKit/USRCAT build: 299,999/299,999 valid, zero failures, 373.93
  seconds single-process. Conformer distribution reconfirmed: 99,999 molecules
  with three and one molecule with two.
- Morgan radius-2/2048-bit baseline: 100,000/100,000 valid, 177.57 seconds.
  Morgan and exact-USRCAT molecule Top-1000 overlap was only five molecules;
  their recall channels are complementary for this out-of-library-type query.
- Real FAISS benchmark used 30,000 training conformers, nlist=512, and a broad
  candidate pool. At nprobe=128, exact USRCAT Top-100/Top-1000 recall was
  1.000/0.989. At nprobe=256, Top-100/Top-1000 were both 1.000 and exact
  Top-10,000 recall in a 100,000-candidate pool was 0.9958. Query latency was
  14--26 ms depending on PQ configuration/nprobe.
- m20, m30, and padded-m32 produced nearly identical broad-candidate recall,
  showing coarse IVF-list coverage—not PQ code length—was limiting at this
  scale. m20 was smallest/fastest and remains in the next sweep.
- RDKit alignment baseline on the exact-USRCAT nearest 1,000 conformers:
  1,000/1,000 valid; pure shape compute 3.33 seconds and jointly optimized
  shape+unweighted-color compute 3.94 seconds. Total wall time was 85.57 seconds
  because the prototype rescanned 3.4 GB MOL2 to locate candidates, directly
  motivating local binary offsets and offset-sorted reads.
- The best RDKit unweighted joint shape+color sum was about 0.58; this low value
  is consistent with QT9 being chemically dissimilar to the macrocycle library
  (maximum Morgan only 0.10). It is not evidence about anchor weighting,
  projected points, Tversky, enrichment, or biological activity.
- Incremental artifact prototype profiling identified RDKit BaseFeatures as the
  dominant cost: on 1,000 conformers, parsing/FeatureFactory/PMI/USRCAT took
  0.75/24.61/0.02/0.14 seconds. A slow half-library shard run was stopped after
  about 30 minutes before atomic promotion; no partial data entered the catalog.
- Feature-template reuse was validated on 100 molecules x 3 conformers: all
  family/atom-membership signatures were identical within each molecule.
- Real multiprocessing benchmark with template reuse, full sanitization,
  features, USRCAT, and PMI: on 3,000 conformers, 1/4/8/16 workers achieved
  95.5/330.2/515.5/647.0 conformers/s with zero failures. On 30,000 conformers,
  8/16/24 workers achieved 413.7/583.8/737.5 conformers/s, again zero failures.
  Twenty-four workers is the current laptop ingestion default, subject to a
  full writer/I/O benchmark; its compute-only extrapolation is about 6.8 minutes
  for 299,999 conformers.
- The 24-worker atomic writer then completed both real shards: 149,999
  conformers in 190.31 seconds and 150,000 in 204.81 seconds (395.12 seconds
  total, about 759 conformers/s including binary writes and checksums). It wrote
  17,174,338 heavy atoms and 12,801,478 base pharmacophore features.
- Both shards had continuous global IDs (0--299,998), exact terminal offsets and
  byte lengths, valid hashes, and bit-identical USRCAT versus the independent
  single-process build (maximum absolute delta 0).
- Quantization audit: 100 random original conformers had maximum heavy-atom
  coordinate error 0.0050004 angstrom. On 100 direct FeatureFactory comparisons,
  stored feature types were exactly identical and maximum center error was
  0.0050004 angstrom, consistent with the 0.01-angstrom grid.
- Warm local mmap reads of complete coordinate+feature records took 9.1 ms for
  1,000 offset-sorted QT9 candidates and 90.5 ms for 10,000; random order took
  13.0 and 110.5 ms. Checksums and record counts agreed. Cold-cache behavior
  remains to be measured separately on the future larger local artifact set.
- The append-only IVF-PQ build trained once on 30,000 records (0.741 s), added
  shard 1 (149,999 conformers) in 0.312 s, persisted/reloaded it, then added
  shard 2 (150,000 conformers) in 0.302 s. The final m20 index contains 299,999
  stable int64 IDs and occupies 8,588,568 bytes; direct-map storage is disabled.
- Incremental-query validation used stage-local exhaustive USRCAT truth. At
  nprobe=256, both one- and two-shard indices recovered 100% of exact top-100
  and top-1,000. Top-10,000 recall was 0.9914 before append and 0.9969 after
  append; QT9 query latency was 9.2 and 13.9 ms. Boundary conformer IDs 0 and
  149,999 self-retrieved after append, confirming stable old/new ID semantics.
- These are warm laptop results and validate append mechanics, not billion-scale
  latency. A later independent MOL2 shard remains required for end-to-end
  production acceptance, but is not needed to continue anchor-scoring work.
- Exploratory direct-contact extraction on 8BJU/QT9 produced three distance-
  supported polar interaction records: ligand HBD atom 30 to CYS379 O at 2.876
  angstrom, ligand HBA atom 29 to CYS379 N at 2.983 angstrom, and ligand HBA
  atom 39 to ASN376 ND2 at 2.838 angstrom. The current mmCIF has no explicit
  hydrogens, so angle is deliberately recorded as null/not-evaluated rather
  than fabricated. These records are provisional ranking evidence; directional
  projection and chemically explicit hydrogen-angle validation remain next.
- Hydrogen preparation is scoped to the query ligand and residues having any
  atom within 5 angstrom. Reduce/Reduce2 is responsible for adjustable hydrogen
  groups and Asn/Gln/His orientation; its output remains ranking evidence only.
- Deterministic 960-point atom-level Shrake--Rupley SASA (1.4-angstrom probe)
  records isolated SASA, complex SASA, and relative burial. QT9's three contacts
  have relative burials 1.000, 1.000, and 0.990 and complex SASA 0.00, 0.00,
  and 0.33 square angstrom. Distance and burial agree; angle remains null
  because Reduce is not installed locally.
- Linux workstation execution is now implemented as a portable handoff: a
  dedicated environment update adds Reduce to `aidd-workstation`; CLI export
  produced a real QT9 pocket containing 25 complete protein residues; and the
  runner preserves the hydrogenated structure, stderr/flip report, command,
  and hashes. Runtime chemical validation remains pending on Linux.

## E022 — Artifact-backed Gaussian batch reranking

Protocol: `experiments/E022-gaussian-artifact-batch-protocol.md`.

Result classification: confirmatory offline orchestration validation; real
QT9 ranking/runtime/enrichment pending.

- Stable global IDs resolve to mmap-backed artifact records without MOL2
  rescans; coordinates use the locked `origin + int16/100` reconstruction.
- Query packages and result packages are hash-bearing NPZ/JSON pairs.
- Centroid/PCA fallback and typed pair seeds recovered locked synthetic rigid
  geometry.
- Every input candidate survived in original order. Shape-only, unweighted
  atom-centered joint, and anchor-weighted atom-centered joint each retained an
  independent best transform and raw overlap primitives.
- Stored scores reproduced at their stored poses and repeated runs were
  numerically deterministic. Full suite: 90 passed.
- Artifact schema v1 cannot support projected color or torsion refinement;
  those remain a versioned companion-artifact task.

## E023 — Staged parallel and resumable Gaussian reranking

Protocol: `experiments/E023-staged-parallel-resumable-gaussian-protocol.md`.

Result classification: confirmatory offline orchestration validation; real
workstation scaling pending.

- The observed 13.2-fold PCA/pair runtime gap motivated a two-stage compute
  policy rather than full-library pair-seed scoring.
- Coarse scoring preserves all input IDs; refine scoring receives the stable
  union of three separately ranked Top-N lists.
- Rigid-invariant self overlaps are cached, candidate chunks run in worker
  processes, and completion order cannot change merged order.
- Atomic, hash-validated chunk manifests provide power-loss recovery. Synthetic
  interruption and one-vs-two-worker checks passed. Full suite: 92 passed.

## E024 - Billion-scale tiered retrieval and Gaussian execution

Protocol: `experiments/E024-billion-scale-tiered-search-protocol.md`.

Result classification: confirmatory offline architecture validation; physical
10M, 100M, and 1B performance validation pending.

- Added exact heap merging of memory-mapped shard-local ranked prefixes.
- Added fixed budgets for baseline, strict, balanced, and loose channels;
  baseline IDs are preserved and duplicate evidence flags are accumulated.
- Added resumable slim coarse Gaussian chunks with a 20-byte logical row and
  exact streaming per-objective Top-N selection.
- Detailed pair refinement remains a separate retained result with a smaller
  refine chunk size so worker concurrency does not collapse on short selections.

## E025 - Multi-cocrystal molecule aggregation before docking

Protocol: `experiments/E025-multi-cocrystal-molecule-aggregation-protocol.md`.

Result classification: confirmatory offline orchestration validation; real
multi-cocrystal workstation inputs pending.

- Each co-crystal query remains an immutable, independently hashed result.
- Conformers collapse to molecules separately for every objective; the primary
  objective also retains configurable Top-M alternative conformers.
- Queries are unioned only inside an explicit site, using within-query ranks
  and RRF rather than incomparable raw-score averages.
- Per-query protected lanes preserve query-specific chemical space before a
  site-consensus lane fills the docking budget.
- Every admitted site/molecule expands into supporting query/receptor-specific
  docking tasks with the selected conformer and transform.

## E026 - Real 8BJU/1X8B dual-cocrystal validation

Protocol: `experiments/E026-real-8bju-1x8b-multi-cocrystal-validation-protocol.md`.

Result classification: confirmatory real Linux workstation validation accepted.

- Query 1 reuses the accepted 8BJU/QT9 detailed result; query 2 is the WEE1A
  1X8B/824 co-crystal at auth chain A/residue 901.
- Both queries are assigned to `WEE1-ATP-site`, while their receptor IDs remain
  distinct so admitted molecules expand into the correct docking hypotheses.
- The runner downloads official coordinate/dictionary inputs, verifies the
  deposited ligand identity, extracts query anchors, reuses the existing
  artifact/pharmacophore/FAISS indices, performs staged Gaussian scoring, and
  invokes E025 aggregation.
- Upstream NPZ products are reused after acceptance rather than regenerated,
  preserving hashes and Gaussian chunk recovery across interruption.
- The first workstation attempt stopped before retrieval: a previously exported
  generic `QUERY_DIR` redirected E026 to the 8BJU manifest, and conversion of
  CCD 824 discarded its explicit Kekule orders. The runner now uses E026-only
  override names and validates the query ID; CCD conversion preserves supplied
  SING/DOUB orders before RDKit aromaticity perception. This was corrected
  before the accepted rerun below.
- The corrected run was accepted. It produced 12,834 union molecules: 955
  shared, 5,991 unique to 8BJU and 5,888 unique to 1X8B. The all-query overlap
  is 7.44% of the union, supporting complementarity rather than redundancy.
- 1X8B contributed five supported anchors, ten invariant pairs, retained all
  299,999 conformers at broad retrieval and refined 8,073 conformers.
- The 5,000 admitted molecules expanded into 5,834 receptor-specific tasks.
  With exactly two queries, this implies 834 shared admitted molecules and
  4,166 single-query admitted molecules.

## E027 - Dual-cocrystal overlap and docking-queue analysis

Protocol: `experiments/E027-dual-cocrystal-overlap-analysis-protocol.md`.

Result classification: implementation validated offline; real E026 report run
pending on the workstation.

- Added hash-verified recomputation of query coverage, union/intersection,
  pairwise Top-K overlap, shared-rank Spearman correlation, objective-specific
  conformer agreement and docking-task multiplicity.
- Emits both `analysis.json` and a human-readable `analysis.md`; no library
  rescoring or molecule reconstruction is required.
- The real workstation report was accepted: Top-100 overlap was zero,
  Top-500/1,000/5,000 intersections were 6/28/505, and the shared-rank
  Spearman correlation across 955 molecules was -0.0371. The 5,000 admitted
  molecules remained 4,166 single-query plus 834 dual-query entries, producing
  exactly 5,834 receptor-specific tasks.

## E028 - Pre-docking aligned-pose and pocket-exclusion QC

Protocol: `experiments/E028-predocking-pocket-qc-protocol.md`.

Result classification: confirmatory offline implementation validation; real
8BJU/1X8B workstation geometry report pending.

- Added deterministic per-query Top-N selection from the immutable E026 docking
  tasks and stable-ID resolution from the one-time conformer artifact catalog.
- Stored 4x4 candidate-to-query transforms are applied without optimization.
  Query coverage, centroid displacement and coordinate-only protein proximity
  are emitted as diagnostic annotations and never as retrieval filters.
- Added generic-element point-cloud PDB exports and a PyMOL review script. The
  files explicitly state that artifact v1 lacks chemistry and that they are not
  valid docking inputs.
- Input/output hashes, receptor identity, query/receptor frame assignment,
  affine-transform validation and repeated-output determinism are enforced.
  The complete offline suite passes 108 tests.

## E029 - One-time chemical/directional companion and fast refinement

Protocol: `experiments/E029-chemical-directional-flexible-companion-protocol.md`.

Result classification: confirmatory offline schema and mathematical-reference
validation; real 299,999-conformer workstation build and timed batch integration
pending.

- Added an append-only per-shard companion keyed by artifact-v1 global IDs. A
  single source pass validates every registry record hash and stores heavy-atom
  identity/charge/aromaticity/chirality, bonds/stereo, feature memberships,
  signed polar directions, axial aromatic normals and bounded terminal torsions.
- Completed shards are reused and relocated source paths are explicit. Query
  execution mmaps the companion instead of rescanning source MOL2.
- Gaussian query packages now retain optional directional arrays; observed
  anchor projection points override geometry-only directions when available.
- Added validated directional Gaussian, element-aware soft vdW exclusion and
  deterministic two-terminal-torsion beam-refinement primitives. E028 consumes
  an optional companion to emit chemical SDF and vdW diagnostics.
- Five new focused checks plus extended query/E028 checks pass; the complete
  dependency-light suite is 115 passed. RDKit production construction remains
  unclaimed until the Linux workstation run.

### 2026-09-08 operational handoff

- The real E029 confirmation remains pending because the Linux path of the two
  original MOL2 shards has not yet been recovered. The E019 artifact shard
  directories intentionally contain no raw MOL2.
- Required sources are `split_0001.mol2` and `split_0002.mol2`; accept them only
  after matching SHA-256 values recorded in
  `to_human/resume-checkpoint-2026-09-08-e029.md`.
- Use explicit `--source split_0001=...` and `--source split_0002=...`
  relocation arguments. This preserves registry source identity while reading
  the physical bytes from university storage or local NVMe.
- Confirmatory measurements remain wall time, peak RSS, worker utilization,
  exact conformer/global-ID counts, source/output hashes, restart behavior, and
  absence of unpromoted partial shards. Directional/torsion batch ranking is a
  subsequent experiment, not part of this source-recovery step.

### 2026-09-09 workstation partial-recovery diagnosis

Result classification: exploratory operational failure analysis; the real E029
build remains confirmatory and pending.

- The two source MOL2 shards are now present under
  `/mnt/local/hand/yuzhang/aidd/mc_data`.
- The workstation build stopped because
  `chemical-companion-v1/.split_0001.partial` already existed. The companion
  format does not implement shard-internal resume, so this directory cannot be
  promoted or appended to without rebuilding the shard.
- Code inspection found that a missing or omitted relocated source path could
  create an empty partial before failing. Preflight now verifies source
  existence, whole-file artifact SHA-256, and registry/artifact identity before
  creating the partial.
- Two regression checks cover missing and hash-mismatched sources without
  partial creation. The complete dependency-light suite passes 117 tests.
- The workstation recovery protocol preserves the old partial outside the
  output root, supplies both explicit source overrides, and captures GNU time
  metrics and a build log. No production build result is claimed yet.

### 2026-09-09 workstation-native identity correction

Result classification: confirmatory implementation correction; real build
still pending.

- The default workstation registry contained zero conformers for the E019
  library. This established that copying the 221 MB Windows test/build registry
  would create an unnecessary and incorrect production dependency.
- E029 now derives stable identity from the immutable artifact-v1 meta and ID
  arrays in source-record order. It requires the relocated Linux MOL2 to match
  the exact whole-file SHA-256 stored in the shard and retains row-level heavy
  atom and feature shape checks during construction.
- Registry validation remains optional for forensic redundancy only. The
  workstation runner no longer passes the active runtime registry and defaults
  to `/mnt/local/hand/yuzhang/aidd/mc_data` for both source shards.
- A new identity-order regression check passes; the full suite is 118 passed.

### 2026-09-09 registry-free completion-count regression

Result classification: exploratory implementation failure; real build remains
pending.

- The first registry-free workstation execution reached the final shard count
  check but raised `UnboundLocalError` because that check still referenced the
  optional registry variable `rows`.
- Replaced the residual reference with the authoritative artifact-v1 row count.
- Added a full registry-free shard regression that exercises binary writing,
  count acceptance, manifest creation, and atomic partial promotion without
  RDKit or SQLite. E029 focused tests pass 5/5; the full suite passes 119/119.
- The failed workstation partial is preserved as diagnostic evidence and must
  be moved outside the output root before the corrected rebuild.

### 2026-09-09 real E029 acceptance and chemical E028 execution

Result classification: confirmatory real workstation validation.

- E029 completed 149,999 and 150,000 conformers in 226.737 and 244.290 seconds:
  299,999 conformers in 471.027 summed shard-seconds, or 636.9 conformers/s.
  Global IDs were contiguous over `[0, 299999)`, four boundary records resolved,
  every manifest file size and SHA-256 matched, and no partial remained.
- Total companion content is 17,174,338 heavy atoms, 18,352,936 bonds,
  12,801,478 features, 18,082,517 feature memberships, 2,334,097 bounded
  terminal torsions, and 14,858,132 torsion moving-atom memberships.
- Identical-command reuse completed in 1.49 seconds with 43,008 KB peak RSS,
  zero filesystem input, and unchanged catalog plus shard manifest hashes. The
  first build's peak RSS and aggregate CPU utilization were not captured.
- Chemistry-aware E028 completed 200 tasks: 100 each for 1X8B/824 and 8BJU/QT9.
  It emitted and hash-validated 200 SDFs, preserved retrieval/admission order,
  and reported chemical topology and vdW annotations as available.
- All 200 poses had points within 2.0 angstrom of receptor atoms and 195 had
  points within 1.5 angstrom. Median candidate/query centroid displacement was
  1.001 angstrom. These are QC annotations, not grounds to delete candidates;
  the prevalence requires distribution analysis, visual review, and local pose
  relaxation during docking/refinement.
- Per-query vdW analysis showed substantially worse rigid-overlay collisions
  for 1X8B/824 than 8BJU/QT9. Median soft penalty/clashing fraction/severe atoms
  were 50.197/0.582/25.5 for 1X8B and 17.337/0.389/15.0 for 8BJU. Soft-penalty
  ranges were 13.942--83.300 and 6.900--49.025, respectively. The highest
  collision poses reached implausible 0.432 and 0.212 angstrom minimum center
  distances. LOW/MEDIAN/HIGH representatives are locked for visual review;
  this analysis remains annotation-only.

## E030 - Fixed-budget chemical pose refinement integration

Protocol: `experiments/E030-fixed-budget-chemical-pose-refinement-protocol.md`.

Status: protocol locked before implementation.

- Human review confirmed correct pocket/frame placement for both queries.
- 8BJU conflicts were mainly terminal with no observed macrocycle penetration.
- 1X8B LOW was terminal dominated; MEDIAN/HIGH cores were too close, supporting
  bounded rigid micro-relaxation in addition to terminal torsions.
- Numeric penalty and visual MEDIAN/HIGH preference did not perfectly agree,
  so named reversible pose variants replace any vdW-only final ranking.
- Implemented the locked 13 rigid micro-seeds: identity, six axis translations,
  and six centroid rotations. Tests confirm deterministic order, exact identity,
  centroid preservation for rotations, and invariant intramolecular distances.
  Batch variant integration remains pending; no accepted ranking changed.
