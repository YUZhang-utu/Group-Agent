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

## 2026-08-28

- Recreated the `aidd-workstation` environment on university Linux storage and
  installed the repository as an editable Python package. Fixed the deployment
  instructions so package installation is explicit after Conda creation.
- Made the optional-RDKit dependency test deterministic across core and
  workstation environments. The workstation suite reached 27 passing tests.
- Restored the approved mandatory Task boundary in the public CLI. Added
  `create-task`, `list-tasks`, `activate-task`, `active-task`, and
  `deactivate-task`; new Campaigns now require matching active Project and Task
  contexts. An end-to-end CLI regression test raised the suite to 28 tests.
- Created the first real university deployment context: user
  `USR-667167B3141D`, Project `PRJ-600CDDACB90A`, Task
  `TSK-A0C50CA9CF19`, and Campaign `CAM-C0B597CE7ADE`.
- Bound the Campaign to human Wee1-like protein kinase (`WEE1`, UniProt
  `P30291`) and executed the first real RCSB search for X-ray structures at no
  worse than 3.0 angstrom resolution.
- The first search exposed a real RCSB GraphQL nullability case for structures
  without nonpolymer entities. Normalized null optional lists, added a
  regression test, and confirmed 29 passing tests.
- The repeated WEE1 search completed and registered 25 PDB candidates. No PDB
  has yet been frozen as the selected receptor.
- Current continuation point: implement or use a terminal-driven PyMOL review
  command to download and visually inspect individual Campaign candidates. Then
  define the ATP-pocket residues, compare multiple structures, and record a
  human-reviewed selection rationale.
- Confirmed on the university environment that a selected PDB candidate can be
  downloaded and opened successfully in PyMOL from the terminal. This completes
  the single-structure GUI smoke test; no final PDB selection was recorded.
- Clarified the next architecture boundary: the existing CLI is the validated
  execution and provenance layer, while Codex natural-language orchestration is
  not yet integrated. The next implementation should map prompts to active
  User/Project/Task/Campaign context and audited tool calls.

## 2026-08-31

- Kept AIDD macrocycle-agent work isolated from the unrelated KRAS
  necessity/enhancement research in the workspace root.
- Completed E009 Campaign-aware multi-structure review integration. Users can
  create a comparison set by PDB ID, inspect resolved members and ligands, and
  generate a Project-scoped multi-structure PyMOL script.
- Added explicit rejection of foreign Campaign candidates, invalid references,
  malformed chains, and absent downloaded mmCIF files. Comparison operations
  preserve `structures_review`; final receptor freezing remains a distinct
  human-reviewed `select-pdb` action.
- Full suite: 31 tests passed.
- Completed E010 after observing that RCSB nonpolymer entities caused GOL, CL,
  NA, EDO, PO4, and MG to appear as ligands. Added conservative component
  classification without discarding provenance, exposed exclusions through
  Campaign status, and made repeated PDB searches refresh existing candidate
  metadata. Full suite: 32 tests passed.
- Completed E011, establishing AI as the prompt-driven recommendation layer
  over the deterministic Project/Campaign execution system. Added immutable,
  provider-neutral AI requests and structured recommendations with evidence
  citation, confidence, uncertainty, model provenance, privacy classification,
  and human-review requirements. Experimental raw data and measurements are
  rejected from model context. Full suite: 35 tests passed.
- Completed E012 predicted target-structure deployment adapters. Added verified
  AlphaFold DB acquisition and unified AlphaFold 2/3, Boltz-2, and Chai-1 input
  and command generation. No GPU inference was claimed locally; real WEE1 runs
  remain a university integration step. Full suite: 38 tests passed.

## 2026-09-01

- Validated the independently deployed AlphaFold 3 Apptainer image on the
  university RTX 5090 workstation with official parameters and the complete
  genetic/template databases.
- Completed ubiquitin and WEE1 kinase-domain predictions. WEE1/PDB PyMOL review
  reported approximately 1.10 Å C-alpha RMSD with all aligned residues and
  0.37 Å after outlier rejection. This supports structural consistency but is
  not an automatic receptor-selection decision.
- Implemented E013 prediction-result registration. Completed AF3 directories
  can now be inspected into hashed manifests, registered beside experimental
  PDB candidates, and selected through a common immutable receptor lock.
- Full offline suite: 41 tests passed. Production WEE1 import remains the next
  workstation action.
- The production WEE1 AlphaFold 3 result was imported into the active Campaign
  after the workstation updated to commit `05e6588`; the user confirmed the
  updated workstation suite passes 41 tests.
- Implemented E014 unified receptor ensembles. All Campaign PDB and predicted
  candidates can be reviewed together with an explicit reference, filtered
  ligand selections, 6 Å ligand pockets, and predicted-construct residue-number
  offsets. Offline suite: 42 tests passed; 8BJU/all-WEE1 PyMOL validation is next.
- Workstation visualization exposed 9TG7 as an incompatible ensemble member.
  RCSB identifies chain A as beta-TrCP and chain B as only a 12-residue WEE1
  degron. Added audited `PDB=reason` ensemble exclusions rather than deleting
  the candidate or silently forcing an alignment.
- Completed E015: provider-neutral LLM receptor-eligibility requests and strict
  response validation now require every candidate, supplied evidence citations,
  and human control. Recommendations do not mutate Campaign state. Full offline
  suite: 44 tests passed.
