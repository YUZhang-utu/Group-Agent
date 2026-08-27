# Research log

## 2026-08-26

- Bootstrapped a separate local-first AIDD agent project to avoid modifying the
  existing WEE1 and lab-os-v4 work.
- Chose batch MOL2 ingestion as the first executable milestone.
- Observed existing MOL2 samples with multiple molecule blocks per file, so the
  importer is explicitly record-oriented rather than file-oriented.
- Locked experiment E001 before validation.
- E001 confirmed: four tests passed and an existing five-record MOL2 file passed
  dry-run validation with no invalid records.
- Deferred the full library import until its exact source directory/file is
  selected; no user structure files were copied or modified.
- Completed dry-run validation of `macrocycles-v1`: 2 source files, 299,999
  records, 0 invalid records, and no issues.
- Completed the production SQLite import as library `LIB-AFA68EE6888C`:
  100,000 molecules, 299,999 conformers, and 34,556,550 total atoms.
- Confirmed grouping distribution: 99,999 molecules have three conformers. One
  molecule (`MOL-3E22C7F32643`, `c--A-Wnme-PdFnme-dW-dLnme-c`) has only conf1
  and conf2; conf0 is absent from the source library. This is a recorded source
  exception, not an importer failure.
- Next continuation point: implement the Campaign state machine, RCSB/PDB
  structure search, human structure selection, and pocket confirmation.

## 2026-08-27

- Locked E002 before implementation. A Git protocol commit was unavailable
  because the project is untracked and the sandbox exposes `.git` read-only.
- Implemented Campaign, target, PDB-candidate, and append-only decision-event
  records with enforced state progression and a frozen selection checkpoint.
- Implemented RCSB search/detail retrieval by UniProt, experimental method, and
  maximum resolution, preserving query and candidate metadata locally.
- Added five Campaign/PDB CLI commands; the full project suite passes 7/7.
- Migrated the production SQLite registry additively. No Campaign was created
  without a target assumption; E003 starts by defining the biological target
  and desired binding/conformational state.
- Completed E003 multi-user task context. Shared models, docking, MD, libraries,
  and caches remain platform resources, while users, projects, tasks, Campaigns,
  decisions, and generated workspaces are isolated.
- Enforced active-task matching for Campaign creation and all Campaign/PDB
  operations. Added user/project/task lifecycle CLI commands and English-only
  task workspace templates.
- Ten tests pass. The additive production migration preserved exactly 100,000
  molecules and 299,999 conformers; no real user or task was invented.
- Completed E004 by simplifying the public hierarchy to user -> active Project
  -> optional Campaign -> Run. Legacy Task storage remains only for migration
  compatibility and is no longer exposed through the CLI.
- Added fixed-root Project workspaces, generic Run/Artifact provenance,
  input/output lineage, metrics, execution metadata, compute keys, scientific
  decisions, and verified result-manifest ingestion.
- Thirteen tests pass. Production migration again preserved exactly 100,000
  molecules and 299,999 conformers. No real storage root or user Project was
  created without an explicit path and identity choice.
- Completed E005 cross-platform Slurm bridge scaffolding. Added Cluster
  Profiles, user cluster identities, OpenSSH and PuTTY command transports,
  validated Slurm Bundle generation, Bundle integrity checks, and content-bound
  submission receipt import.
- Seventeen offline tests pass across platform commands, registered-root path
  safety, Slurm validation, tamper detection, and receipt matching. Production
  migration preserved the molecular registry and created no fictitious cluster
  records or external jobs.
- Began E006 multi-structure target review and ligand similarity. Added named
  multi-candidate comparison sets, pocket Kabsch alignment with residue-level
  deviations and missing coverage, Project-scoped PyMOL scripts, and optional
  RDKit Morgan/USRCAT adapters.
- Twenty-one tests pass. The current laptop lacks RDKit, Gemmi, BioPython, and
  PyMOL, so persistent similarity-index construction, mmCIF parsing, and GUI
  launch remain deployment validations rather than claimed results.
- Locked E007 before implementation and added separate reproducible core and
  workstation Conda environments for Python 3.11, RDKit, BioPython, Gemmi, and
  PyMOL, plus a machine-readable environment diagnostic.
- Added validated project-scoped RCSB mmCIF acquisition with SHA-256 metadata.
  Added AlphaFold 3 input generation and model-profile validation while keeping
  model parameters, databases, credentials, and predictions outside Git.
- Added disk-backed Morgan fingerprint and USRCAT descriptor indices with
  version validation and deterministic ranking. Twenty-seven tests pass; the
  dependency-backed RDKit and AlphaFold runtime checks remain for the desktop.
- Audited source-control scope: local libraries, SQLite registries, structures,
  indices, model weights, prediction outputs, secrets, and workspaces are
  ignored. The uploadable code/documentation set is approximately 181 KB.
