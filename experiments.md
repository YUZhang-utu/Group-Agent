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

Status: active — confirmatory integration

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

## E017 — Hierarchical 2D/3D macrocycle-library search

Status: planned — confirmatory integration

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
