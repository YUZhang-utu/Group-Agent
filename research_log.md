# Research log

## 2026-09-02

- Resumed only the Group-Agent AIDD macrocycle project from the 2026-09-01
  checkpoint; no KRAS necessity/enhancement files were read or modified.
- Locked E018 before implementation.
- Added retained-instance mmCIF enumeration with model, author chain/residue,
  insertion code, and altloc identity.
- Added validated RCSB CCD caching and CCD-bond-based RDKit reconstruction of
  crystal-coordinate ligands. No distance-based bonds or generated conformers
  are permitted; missing CCD heavy atoms are visible failures.
- Added versioned parent standardization, SDF/USRCAT/manifest output, and CLI
  integration for direct Campaign ligand registration.
- Added a source-streaming production builder for molecule-level Morgan and
  conformer-level USRCAT indices, the conformer mapping, and checksummed
  manifests for registered MOL2 libraries.
- Desktop validation: 53 tests passed. RDKit and Gemmi are absent from the
  desktop environment, so real WEE1 chemistry and the 100,000-molecule build
  remain explicitly pending on the university workstation.
- Accepted a revised billion-scale search design and locked E019 before its
  implementation. The revision makes conformer coverage a target-independent
  calibration on Platinum, requires query crystallographic/protonation QC,
  uses projected protein-interaction anchors, and evaluates ColorTversky.
- Corrected unsafe index pruning: key admission is now governed only by global
  document frequency and type rules, never a per-conformer rare-key cap.
- Decoupled distance tolerance from bin width, added degenerate-triangle and
  Gaussian cutoff rules, and made brute-force top-N index recall mandatory.
- Deferred native SIMD/CUDA implementation until a one-million-scale prototype
  measures index recall, posting distributions, enrichment ablations, and
  actual stage latency.
- Refined E019 with an exact two-pass 4096-bucket global-df computation. Only
  high-df keys are excluded; all rare keys remain eligible, and admitted-key
  postings are generated in a second pass. Added shared-key survival >=90% as
  a direct admission-rule diagnostic.
- Specified directional HBD/HBA projection geometry, consistent hydrogenation
  and protonation manifests, projected-vs-atom-centered bin regimes, explicit
  q-biased Tversky fields, and renamed the weighted score `anchor_combo` to
  avoid false comparison with ROCS ComboTanimoto.
- Corrected multi-projection weighting: under product-amplitude overlap,
  `W/sqrt(n)` per point preserves functional-feature self-overlap; `W/n` does
  not. Added fixed-N hit rate and the full Morgan-to-anchor baseline ladder.
- Corrected the E019 execution environment: billion-scale df counting and
  index construction run on the local desktop, not a Slurm cluster. The 4096
  buckets are processed by a bounded local worker pool with NVMe-aware I/O and
  per-bucket checkpoints.
- Replaced the E019 triplet/df/bin index route entirely with an in-memory FAISS
  IVF-PQ USRCAT L1 and local-NVMe bounding-box/PMI L1.5. NFS is now read-only
  cold input scanned once per shard; uncompressed coordinates, features,
  metadata, PMI/bounds, and raw USRCAT are emitted together and checkpointed.
- Preserved all supplied conformers and retained raw float32 USRCAT vectors for
  cheap FAISS retraining. Flagged two production constraints for validation:
  60-D PQ cannot directly use m=32, and one-billion-vector FAISS IDs add about
  8 GB beyond PQ codes.
- Split chemistry perception from fast conformer coordinate ingestion because
  unsanitized RDKit records cannot reliably define aromaticity, hybridization,
  formal charges, hydrogens, or pharmacophore types.
- Raised the m30 FAISS budget to 42--46 GB steady-state and 52--56 GB including
  workspaces; required direct_map-off validation and clean-process reload after
  index construction. Standard FAISS int64 IDs remain unless a custom implicit
  layout is separately validated.
- Added m20/m30/padded-m32 PQ ablations and per-dimension USRCAT z-scoring;
  optional polar-block weighting occurs after z-score.
- Added source grouping/identity audit, conformer-specific projection geometry,
  parent-feature count pruning, cold/warm and offset-sorted read benchmarks,
  held-out threshold calibration, and result transforms with top-1000 SDF
  export.
- Ran the first laptop-scaled E019 FAISS smoke tests without touching NFS or the
  full MOL2 library. At one million synthetic 60-D vectors with nlist=1024,
  nprobe=8 and exact top-100 evaluation, candidate recall was 0.9010 (m20),
  0.9464 (m30), and 0.9474 (zero-padded m32). Serialized sizes were 28.32,
  38.32, and 40.34 bytes/vector respectively. The result provisionally favors
  m30 over padding but remains exploratory until repeated on real USRCAT.
- The strict one-million run completed in 16.45 seconds with maximum reported
  process RSS 0.864 GiB. This validates safe laptop-scale testing only and does
  not include chemistry ingestion or Gaussian L2 timing.
- Selected WEE1 8BJU/QT9 as the first real query after validation QC: 1.53
  angstrom, ligand RSCC 0.933, EDIA 0.951, occupancy 1.0, all 42 atoms supported,
  and CCD-authoritative bond orders.
- Built real USRCAT for all 299,999 registered conformers in 373.93 seconds:
  every record succeeded. Built Morgan radius-2/2048 for all 100,000 molecules
  in 177.57 seconds with zero failures. Their Top-1000 overlap was only five.
- Corrected the FAISS metric from approximate rank agreement to exact Top-N
  containment in a broad L1 pool. On QT9, nprobe=128 recovered exact USRCAT
  Top-100/1000 at 1.000/0.989; nprobe=256 reached 1.000/1.000 and recovered
  0.9958 of exact Top-10,000 in 100,000 candidates at 21--26 ms.
- Ran RDKit shape-only and jointly optimized unweighted shape/color alignment
  on 1,000 real conformers. Compute was 3.33/3.94 seconds, but rescanning MOL2
  made wall time 85.57 seconds, validating the need for local offset-addressed
  binary artifacts. Anchor/Tversky and enrichment conclusions remain pending.
- Implemented the first immutable append-only conformer artifact shard schema:
  int16 heavy-atom coordinates, typed int16 pharmacophore centers, binary
  offsets/IDs, bbox/PMI, feature counts, raw USRCAT, hashes, and stable global
  IDs suitable for later FAISS `add_with_ids`. New library material appends a
  shard instead of rewriting old artifacts.
- Aborted the first 149,999-conformer full-feature shard after about 30 minutes
  because it had not reached atomic completion. A 1,000-record profile showed
  24.61 seconds in RDKit FeatureFactory versus 0.75 parse, 0.14 USRCAT, and
  0.02 PMI seconds; generic pharmacophore matching was 96.4% of measured core
  time. The hidden partial shard was never promoted or cataloged.
- Verified on 100 consecutive three-conformer molecules that feature family and
  atom-membership templates were identical across all three conformers. Updated
  the builder to compute the graph feature template once per molecule and
  recompute only feature coordinates for each conformer. Multiprocess scaling
  remains the next ingestion benchmark before another full shard run.
- Completed real bounded multiprocessing tests after template reuse. All tested
  records succeeded. At 3,000 conformers, 1/4/8/16 workers delivered
  95.5/330.2/515.5/647.0 conformers/s. At 30,000 conformers, 8/16/24 workers
  delivered 413.7/583.8/737.5 conformers/s. The current laptop sweet spot is 24
  workers, implying about 6.8 minutes compute time for the 299,999-conformer
  library before writer/I/O overhead is measured.
- Integrated 24-worker processing into the atomic shard writer and completed
  both production-scale local shards: 299,999 conformers in 395.12 seconds,
  including uncompressed writes and hashes. Global IDs, offsets, file lengths,
  and hashes passed; raw USRCAT was bit-identical to the independent build.
- Audited int16 precision against source MOL2: maximum heavy-atom error across
  100 random conformers was 0.0050004 angstrom. Stored feature types matched
  direct RDKit FeatureFactory exactly across 100 records, with the same maximum
  coordinate error.
- Replaced L2 MOL2 rescans with mmap offset reads. Complete coordinates/features
  for 1,000/10,000 QT9 candidates loaded in 9.1/90.5 ms when offset-sorted and
  13.0/110.5 ms in random order under the current warm-cache condition.
- Superseded the earlier `W/sqrt(n)` projection proposal with the user's locked
  normalized-mixture convention `W/n` per alternative direction. Its changed
  self-overlap is accepted explicitly and all score thresholds must therefore
  be calibrated rather than borrowed from ROCS.

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
- Implemented E016 and the E017 orchestration core. Campaign co-crystal ligands
  are now filtered, provenance-rich instances; query selection is immutable and
  human-controlled; LLM recommendations remain advisory. Hierarchical searches
  retain Morgan, best-conformer USRCAT, 2D fallback, index hashes, and query
  snapshots. Full offline suite: 49 tests passed. Production RDKit/index
  validation remains on the university workstation.
- Recorded continuity checkpoint `to_human/resume-checkpoint-2026-09-01.md` at
  commit `2859048`. The next session begins with automatic mmCIF/CCD ligand
  extraction and production Morgan/USRCAT index construction, not additional
  receptor-selection work or the unrelated KRAS necessity/enhancement project.
## 2026-09-02 - E019 append-only FAISS validation

- Built and reloaded an IVF-PQ checkpoint after each immutable real-library
  artifact shard; counts progressed from 149,999 to 299,999 without rewriting
  the first shard or changing its global ID range.
- Confirmatory QT9 validation at nprobe 256 recovered exact stage-local USRCAT
  top-100/top-1,000 fully. Top-10,000 recall was 0.9914 for shard 1 and 0.9969
  after shard 2; warm query latency was 9.2/13.9 ms.
- Boundary IDs 0 and 149,999 self-retrieved from the final index. This confirms
  the append mechanism on existing data; a newly supplied MOL2 batch will later
  serve as an independent end-to-end acceptance test.
- Full dependency-light test suite: 60 passed. The chemistry environment lacks
  pytest, so tests were correctly run with the base environment and `src` on
  PYTHONPATH; chemistry/FAISS validation used the dedicated mol environment.

## 2026-09-02 - E019 anchor recall-safety amendment

- Corrected the distinction between an observed crystal hydrogen bond and an
  unknown biologically necessary interaction. Anchor evidence may change
  ranking but cannot remove an L1/L2 candidate.
- Withdrew zero-weight solvent features and anchor-count hard gates. Default
  soft weights are now 1.5/1.0/0.5 for anchor/ordinary/solvent-exposed sites.
- Added a lossless candidate-score schema with atom-centered, projected,
  unweighted, anchored, and per-anchor outputs plus reversible reranking and
  Spearman/top-N overlap diagnostics.

## 2026-09-02 - E019 reproducible anchor/pose schema

- Added a query-manifest schema that snapshots stable anchor IDs, ligand atom
  indices, types, atom/projected coordinates, hydrogen-bond/burial evidence,
  and materialization weights.
- Replaced the ambiguous single transform with named objective-specific poses.
  Each pose owns its transform and complete scores evaluated at that transform.
- Stored unnormalized per-anchor cross-overlap and query self-overlap arrays;
  the reranking API resolves them to manifest anchor IDs for readable output.
- Verification: 64 dependency-light tests passed.

## 2026-09-02 - QT9 hydrogen-independent anchor evidence

- Implemented deterministic atom-level Shrake--Rupley SASA and stored isolated,
  complex, and relative-burial values in the query anchor snapshot.
- QT9's three provisional contacts all satisfy the observed 2.6--3.5 angstrom
  heavy-atom window and have 0.990--1.000 relative burial. These concordant
  facts strengthen interaction evidence without asserting necessity.
- Reduce is not installed on this Windows laptop. Angle remains null and
  Asn/Gln/His flipping remains pending; no fallback geometry was presented as a
  Reduce result. Tests: 66 passed.

## 2026-09-02 - Linux Reduce workstation handoff

- Added a Linux update specification for the existing `aidd-workstation`
  environment with Bioconda Reduce, without making it a Windows dependency.
- Added deterministic export of the query ligand plus complete protein residues
  contacting it within 5 angstrom. The real 8BJU/QT9 export contains 25 protein
  residues and was written successfully on Windows for transfer/testing.
- Added a safe argument-array Reduce runner and immutable manifest containing
  scope, command, hashes, hydrogenated PDB, and stderr/flip report.
- Offline runner tests and the complete suite pass: 68 tests. Real hydrogen/
  flip/angle validation is pending execution on the Linux workstation.

## 2026-09-02 - Workstation dependency hotfix

- Linux validation exposed an undeclared top-level `psutil` import in the FAISS
  benchmark. Moved it into the RSS-measuring execution path so pure helpers and
  test collection do not require optional monitoring dependencies.
- Declared both `psutil` and `faiss-cpu` in the standard and Linux workstation
  environments and the Python workstation extra. Full suite: 69 passed.

## 2026-09-02 - Reduce output acceptance gate

- Clarified that a null HET dictionary means Reduce default lookup, not absence
  of a dictionary. Added validation of manifest hashes, ligand/protein hydrogen
  count increases, dictionary/connectivity warnings, and flip records.
- Angle calculation is enabled only when those checks pass. The validator writes
  `reduce_validation.json`; all evidence remains query annotation/reranking and
  does not affect retrieval. Full suite: 70 passed.

## 2026-09-02 - First real QT9 anchor snapshot

- Implemented an mmCIF protein-contact extractor that labels direct polar
  contacts without turning them into candidate filters.
- The initial 8BJU/QT9 snapshot contains three contacts at 2.84--2.98 angstrom.
  Hydrogen-bond angle is explicitly null because the structure lacks hydrogens;
  directional validation/projection is still pending and no necessity claim is
  made from these distance-supported observations.
- Query weights and evidence are serialized beside stable anchor IDs in the
  ignored local `query_manifest.json`. Verification increased to 65 tests.

## 2026-09-02 - QT9-specific Reduce dictionary correction

- Linux reported zero ligand hydrogens because the installed default Reduce HET
  dictionary lacks QT9. This is a fatal ligand-angle failure, not an acceptable
  null-dictionary case.
- Added authoritative CCD-to-Reduce dictionary generation. The real QT9 CCD
  produced 76 atom definitions, including 34 hydrogens and 81 bonds.
- Split validation warnings into fatal connectivity/dictionary problems and
  expected pocket-fragment "appear unbonded" warnings. Full suite: 71 passed.

## 2026-09-02 - Reduce pre-existing hydrogen acceptance correction

- Corrected the validator: a valid output needs ligand and protein hydrogens,
  not a positive before-to-after count delta. Reduce may preserve/reorient an
  already hydrogenated ligand without changing its hydrogen count.
- Hydrogen deltas remain diagnostic fields. Added a regression case with one
  pre-existing ligand hydrogen and zero ligand delta; it is accepted when the
  output contains ligand/protein hydrogens and no fatal warning. Tests: 73.

## 2026-09-02 - Reduce angle and projection enrichment

- The Linux QT9 run passed after using the official wwPDB legacy HET entry and
  added all 34 ligand hydrogens. Protein hydrogenation and GLN375 flip evidence
  are therefore suitable for directional annotation.
- Added query-manifest enrichment from the hydrogenated PDB: donor hydrogen,
  D-H...A angle, H...A distance, and observed protein-partner projection point
  are stored per stable anchor ID. Missing hydrogens remain explicitly null.
- This enrichment changes evidence/reranking only and cannot remove candidates.
  Full suite: 74 passed.

## 2026-09-02 - Continuity checkpoint after Linux Reduce validation

- User confirmed the official wwPDB QT9 HET entry produced an accepted Reduce
  run with all 34 ligand hydrogens. This supersedes the failed default/generated
  dictionary attempts while preserving them as useful negative results.
- Recorded the exact implementation/validation boundary and next Linux command
  in `to_human/resume-checkpoint-2026-09-02.md`. Real angle/projection enrichment
  remains pending; code availability is not misreported as completed evidence.

## 2026-09-02 - Reduce explicit dictionary activation hotfix

- The first Linux rerun added 160 protein hydrogens and evaluated a GLN375
  flip, but added zero QT9 hydrogens; stderr proved the configured dictionary
  was not actually loaded through environment lookup.
- The runner now passes `-DB <absolute dictionary path>` explicitly, verifies
  the dictionary exists before execution, and records
  `dictionary_cli_applied` in its manifest. Tests: 72 passed.

## 2026-09-06 - End-to-end DMTA and laboratory-automation roadmap

- Extended the architecture roadmap from retrieval through validated docking,
  independent rescoring, versioned activity/property prediction, Pareto
  nomination, inventory/sample identity, assay execution, QC, and active
  learning. Current E019 scope and completion claims remain unchanged.
- Defined laboratory automation as a controlled execution plane. The LLM is
  advisory and cannot send free-form device commands; physical Runs require a
  typed/versioned protocol, deterministic device compilation, capability and
  safety checks, simulation or dry-run, and explicit human release.
- Added phased gates so sample lineage and a manually validated assay digital
  twin precede hardware control, supervised execution precedes unattended
  operation, and QC-approved prospective data precede closed-loop retraining.
- Recorded the roadmap in `docs/end-to-end-roadmap.md`; this is an architecture
  decision, not an executed laboratory experiment or validation claim.

## 2026-09-06 - Reusable partial 3D pharmacophore retrieval

- Locked E020 before implementation. The hypothesis separates an expensive
  once-per-library shard index from small per-co-crystal query compilation.
- Added immutable pharmacophore feature-type-pair/distance-bin postings built
  directly from existing `meta.bin` and `feats.bin`; source MOL2 is not read.
- Added library-independent hashed query plans and nested loose, balanced and
  strict partial-motif retrieval. Results form a lossless union with external
  FAISS/L1 candidates, preserving the anchor reranking-only safety policy.
- Fixed Windows atomic promotion by releasing mmap handles before renaming the
  completed shard index. Added CLI and operating documentation.
- Offline synthetic validation and the complete suite pass: 78 tests. A real
  full-library QT9 build/recall/runtime run remains pending on Linux.
- Added a single-command workstation validator that can generate a fresh QT9
  FAISS L1 set from the persisted incremental index, build or reuse the local
  pair index, compile the query, execute all tiers, verify nesting/ranges/counts
  and prove every external L1 ID survived. It writes a machine-readable
  acceptance report and retains all intermediate inputs and outputs.

## 2026-09-07 - Gaussian overlay mathematical reference

- Locked E021 before implementation to prevent score definitions from changing
  after observing real QT9 rankings.
- Implemented deterministic weighted Gaussian shape/color primitives, raw
  cross/self preservation, standard Tanimoto, directed query-biased Tversky,
  parent-feature alternative-point weight conservation, candidate-to-query
  homogeneous transforms, weighted Kabsch, and compatible pair seeds with six
  axial rotations.
- Confirmatory synthetic tests recovered known proper rigid geometry, enforced
  same-type color matching and rejected invalid numerical inputs. Nine focused
  tests and the full dependency-light suite pass: 87 tests.
- This completes the mathematical reference kernel only. Real artifact batch
  integration, objective-specific optimization and QT9 enrichment remain the
  next experiment.

## 2026-09-07 - Artifact-backed Gaussian batch reranking

- Locked E022 before implementation. The confirmatory hypothesis required
  stable-ID artifact reads, lossless L1 membership, deterministic rigid seeds,
  and independent transforms for all named objectives.
- Added read-only catalog/shard resolution and exact coordinate dequantization;
  stale cross-platform absolute shard paths fall back to catalog-local names.
- Added a hashed co-crystal Gaussian query package, principal-axis fallback
  seeds, and batch scoring for shape-only, atom-centered unweighted joint, and
  atom-centered anchor-weighted joint objectives.
- Synthetic artifact tests recovered a locked rigid pose, reproduced scores at
  stored transforms, preserved candidate order, and were deterministic across
  reruns. Full dependency-light suite: 90 tests.
- Result is confirmatory for orchestration correctness only. Real QT9 runtime,
  ranking and enrichment remain pending on the workstation. Projected color,
  exclusion volume, and terminal torsion are not claimed by artifact schema v1.

## 2026-09-07 - Staged parallel and resumable Gaussian execution

- User workstation confirmed the E022 real runner on 1,000 candidates. Pair
  mode required 219.70 seconds versus 16.59 seconds for PCA-only; CPU was
  approximately one logical core and RSS approximately 55 MB in both runs.
- Locked E023 before implementation: score all candidates coarsely, refine the
  deterministic union of three objective Top-N sets, and preserve atomic chunks
  across interruption or power loss.
- Added invariant self-overlap caching, process-local mmap/query initialization,
  configurable workers/chunks/progress, stable merge, and run-configuration
  hashing. Added coarse/refine/all CLI stages.
- Offline confirmatory tests showed exact one-vs-two-worker arrays, complete
  coarse ID retention, exact Top-N union, corrupt-chunk repair without rewriting
  valid chunks, configuration mismatch rejection, and no leftover partials.
  Full suite: 92 passed.
- Real 16-worker throughput and interruption recovery remain to be measured;
  no linear speedup is inferred from synthetic tests.

## 2026-09-07 - E023 real full-library staged run

- Workstation completed PCA coarse scoring for 299,999/299,999 conformers in
  40.83 seconds inside the runner using 300 chunks and 16 requested workers.
- The union of three independently ranked Top-5,000 lists contained 7,704
  conformers, substantially below the 15,000 disjoint-list maximum.
- Pair refinement completed 7,704/7,704 in 77.99 seconds inside the runner;
  whole-command wall was 79.67 seconds, CPU 743%, peak RSS 171,356 KiB, with no
  swaps or major page faults.
- Refine had only eight 1,000-record chunks, structurally limiting active
  concurrency to eight workers and explaining the measured CPU utilization.
- Run manifest status is complete and all stages/results have content hashes.
  Next checks are score monotonicity/rank overlap, molecule-level collapse,
  pose review, and a real resume-reuse invocation.

## 2026-09-07 - E024 scalable tiered execution implemented

- Locked E024 before implementation and separated evidence retention from
  expensive-compute admission.
- Implemented deterministic global Top-K merging over memory-mapped ranked
  shard prefixes. Tie order is score descending then stable global ID.
- Implemented a fixed-budget baseline/strict/balanced/loose schedule. The
  baseline prefix is asserted unchanged; later channels only add IDs or flags.
- Implemented sharded slim PCA coarse output containing int64 ID plus three
  float32 objectives, exact streaming Top-N union, and detailed pair refinement.
- Added distinct coarse/refine chunk sizes, atomic hash validation, resume, CLI,
  operating documentation, and synthetic equivalence checks against E023.
- Offline validation is confirmatory for algorithms only. No billion-scale
  latency, recall, or throughput claim is made before the physical scale ladder.

## 2026-09-07 - E025 multi-cocrystal pre-docking aggregation

- Locked the E025 protocol before implementation.
- Added a versioned plan that binds every query to an exact receptor, biological
  site, detailed Gaussian result, optional manifest and content hashes.
- Added objective-specific conformer-to-molecule collapse with deterministic
  ties and retained Top-M primary conformers.
- Added site-scoped cross-query union, support counts, rank percentiles and RRF;
  raw Gaussian scores are never averaged across query ligands.
- Added protected per-query admission lanes followed by a consensus lane, plus
  receptor/query-specific docking task records carrying conformer and transform.
- Offline synthetic validation passes as part of the 102-test suite. Real WEE1
  multi-cocrystal aggregation remains the next confirmatory workstation run.

## 2026-09-07 - E025 real QT9 single-query smoke accepted

- User confirmed the corrected complete validator finished without error on the
  real QT9 E025 output.
- Accepted checks covered result/source hashes, unique query-molecule and
  site-molecule identities, objective-specific 4x4 transforms, Top-M conformer
  order, admission uniqueness/reasons, docking task support, and manifest counts.
- `raw_scores_averaged=false` is an intentional invariant, not a failed check.
- This promotes the conformer-to-molecule and docking-task data path to real
  single-query validated. Cross-query RRF and protected-lane behavior still need
  at least two independent co-crystal search results.

## 2026-09-07 - E026 real 8BJU/1X8B validation prepared

- Locked the protocol before implementation. The independent queries are
  8BJU/QT9 and 1X8B/824 auth A:901, assigned to one WEE1 ATP site but retaining
  separate receptor identities.
- Added a CLI boundary for atomic co-crystal anchor-manifest extraction and an
  end-to-end workstation runner that reuses the accepted QT9 result and all
  library-scale indices.
- The runner retains accepted retrieval/query NPZ files across invocations so
  staged Gaussian chunk hashes remain stable and power-loss resume is real.
- Offline suite passes 103 tests. No 1X8B retrieval, overlap, RRF, or docking-task
  result is claimed until the Linux workstation run completes.

## 2026-09-07 - E026 first workstation attempt failed before scoring

- The log showed reuse of `data/e019_query_8bju/query_manifest.json`, proving an
  exported generic `QUERY_DIR` had overridden the E026 default. The runner now
  ignores generic output-directory variables, accepts only `E026_*` overrides,
  and asserts the persisted query ID is exactly `1X8B:824:A:901`.
- RDKit then raised `KekulizeException` for CCD 824. Root cause was conversion
  code replacing every bond marked aromatic by wwPDB with `BondType.AROMATIC`,
  even though the dictionary supplies an explicit alternating SING/DOUB Kekule
  assignment. Preserving those authoritative orders sanitized successfully on
  the complete official 824 heavy-atom graph and produced the expected neutral
  fused `[nH]` structure.
- This is an infrastructure-negative result, not a retrieval or enrichment
  result. No library-scale stage ran. The hotfix suite passes 104 tests and a
  clean workstation rerun is pending.

## 2026-09-07 - E026 real dual-cocrystal validation accepted

- The corrected 8BJU/QT9 plus 1X8B/824 run completed with `accepted=true` and
  preserved the no-cross-query-score-averaging invariant.
- The molecule union contains 12,834 entries: 955 shared, 5,991 8BJU-only and
  5,888 1X8B-only. Shared molecules are 7.44% of the union and approximately
  13.75%/13.95% of each query's retrieved molecule set, so the second query
  adds substantial nonredundant chemical space.
- 1X8B used five supported anchors and ten anchor pairs; broad retrieval kept
  299,999 conformers and detailed refinement selected 8,073.
- 5,000 admitted molecules generated 5,834 receptor-specific docking tasks.
  Under the two-query invariant, 834 admitted molecules carry both receptor
  hypotheses. This is confirmatory for orchestration and complementarity, not
  for binding or docking accuracy.
- Locked E027 and implemented a hash-verified Top-K/rank/pose-agreement/queue
  analysis report. Real report generation remains pending on the workstation.

## 2026-09-08 - E027 real overlap analysis accepted

- The real E026 analysis reproduced 12,834 union and 955 shared molecules.
- Pairwise intersections were 0/6/28/505 at Top-100/500/1,000/5,000.
- Shared-molecule rank Spearman was -0.0371, supporting independent protected
  query lanes rather than intersection filtering or raw-score averaging.
- The queue identity closed exactly: 4,166 single-query plus 834 dual-query
  admissions generated 5,834 receptor-specific tasks.

## 2026-09-08 - E028 pre-docking geometry QC implemented offline

- Locked E028 before implementation.
- Added stable-ID candidate materialization, stored-transform application,
  query-frame coverage, receptor-chain protein distance diagnostics, generic
  point-cloud PDB export, PyMOL review generation and hash-bearing manifests.
- The operation reads the immutable artifact catalog and accepted aggregation;
  it does not rescan source MOL2, rescore the library or mutate admission.
- Artifact v1 lacks elements/bonds/topology, so exports are explicitly geometry
  only and prohibited as docking inputs. Chemical SDF export remains assigned
  to a once-built artifact v2 companion.
- Two focused tests and the complete 108-test suite pass. Real 8BJU/1X8B E028
  execution remains pending on the Linux workstation.

## 2026-09-08 - E029 chemical companion foundation implemented offline

- Locked E029 before implementation and retained artifact v1 as the immutable
  recall source.
- Added a multiprocessing, ordered, per-shard companion builder. It performs one
  sequential source-MOL2 pass, checks every conformer against its registry hash
  and artifact-v1 identity, and writes mmap-ready chemistry, topology,
  directional-feature and terminal-torsion arrays.
- Added explicit source relocation overrides so moving a source file does not
  change registered identity. Completed shards are reused; partial shards are
  never deleted automatically.
- Added signed polar/axial aromatic directional Gaussian overlap, element-aware
  vdW penetration, and deterministic bounded terminal-torsion beam refinement.
  Query packages now persist directions and use observed protein projection
  points for anchors when available.
- E028 optionally consumes the companion to produce chemical SDF files and vdW
  diagnostics while retaining the original generic point clouds and immutable
  admission order.
- The full local suite passes 115 tests. Real RDKit construction, throughput,
  binary identity audit and query-time batch timing remain pending on Linux.

## 2026-09-08 - E028/E029 continuity checkpoint

- Pushed E028/E029 protocol and implementation commits through `b401393` to
  `origin/main`; the local and remote `main` branches matched at that commit.
- Recorded E026/E027 real dual-cocrystal acceptance as the current retrieval
  baseline: 12,834 union molecules, 955 shared, 5,000 admitted, and 5,834
  receptor-specific docking tasks. The low overlap continues to require a
  protected union of independently searched co-crystals.
- E028 is implemented offline for aligned-pose point-cloud export and
  pocket-exclusion diagnostics. E029 is implemented offline as a one-time,
  global-ID-keyed chemistry/topology/direction/torsion companion. The complete
  dependency-light suite remains 115 passed; real 299,999-conformer E029
  construction is not yet claimed.
- Clarified the source-data boundary. E019 shard directories contain processed
  binary coordinates/features/USRCAT/IDs, not original MOL2 records and not
  enough authoritative topology for E029. The registered Windows sources are
  `6aa/split_0001.mol2` and `6aa/split_0002.mol2`, with SHA-256 values
  `bc729d337ba7065ebeccdbd7aa05f2985d4ecc778dd8b7949e1b57ada49a3b2b`
  and `6458081fa50eb2676d818497403084775fa13fbc87cbcfbe6c7318f40141d1e7`.
- The user remembers uploading both sources to university storage, but the
  current path has not been recovered. If necessary, transfer them to
  `/mnt/local/hand/yuzhang/aidd/source/6aa/`. E029 already accepts repeatable
  `--source SHARD_NAME=/relocated/source.mol2` overrides, so relocation does not
  require another code change or GitHub commit.
- Locked the next confirmatory action: verify both source hashes, build the
  companion on the Linux workstation, record time/RSS/reuse/counts/hashes, then
  run E028 with chemical SDF and vdW diagnostics. Batch integration of
  directional and terminal-torsion refinement follows only after this real
  build is accepted.
- Based on the measured 299,999-conformer writer, a 300-million-conformer
  artifact build linearly extrapolates to about 4.6 days. The practical
  already-3D 100-million-molecule budget is provisionally 1--2 weeks, while the
  fixed-budget per-query target is 2--10 minutes. Both remain estimates pending
  the prescribed 10M/100M physical scale gates.
- AlphaFold 3 ligand-complex exploration performed in the same conversation is
  an independent workstation task and is intentionally excluded from the
  AIDD-agent experiment and necessity/enhancement evidence chain.

## 2026-09-09 - E029 stale-partial failure diagnosed and preflight hardened

- The user located both original MOL2 shards at
  `/mnt/local/hand/yuzhang/aidd/mc_data` and reported that the companion build
  was blocked by `.split_0001.partial`.
- Determined that companion shards are atomic all-or-nothing builds; their
  partial directories are diagnostic remnants rather than resumable state.
- Found and fixed a preflight ordering defect: the builder created the partial
  directory before validating the relocated source path. It now checks source
  existence, exact whole-file SHA-256, and registry/artifact identity before
  creating output state.
- Added regression tests for missing and hash-mismatched sources. Focused E029
  checks pass 3/3 and the complete dependency-light suite passes 117/117.
- Wrote a recoverable workstation procedure that inventories and moves the old
  partial instead of deleting it, then rebuilds with explicit `mc_data` source
  overrides and records GNU time output. Real E029 acceptance remains pending.

## 2026-09-09 - Removed Windows registry from E029 workstation execution

- Workstation diagnostics showed that its active 544 KB registry has zero
  conformers for `LIB-AFA68EE6888C`; the matching 299,999-row registry existed
  only as a 221 MB Windows test/build artifact.
- Corrected the deployment model instead of copying or overwriting registries.
  E029 now uses exact relocated MOL2 whole-file hashes and artifact-v1 ordered
  global/conformer/molecule identities, with row-shape validation during build.
- Made `--db` optional redundant validation and removed it from the workstation
  runner. The runner now defaults both sources to the Linux `mc_data` directory.
- Added a stable identity-order test; E029 focused tests pass 4/4 and the full
  dependency-light suite passes 118/118. Real build timing remains pending.

## 2026-09-09 - Fixed registry-free final count regression

- The first workstation-native build exposed one residual reference to the
  optional SQLite `rows` variable at final count acceptance.
- Changed the invariant to compare written rows against `len(v1_meta)`, the
  authoritative count in both registry-backed and registry-free modes.
- Added an end-to-end mocked shard build covering output files, final count,
  manifest generation and atomic promotion. Focused tests pass 5/5 and the full
  suite passes 119/119.

## 2026-09-09 - Accepted real E029 and executed chemistry-aware E028

- Accepted the real 299,999-conformer companion after source/output hashes,
  contiguous global IDs, four boundary reads, and absence of partials passed.
  Shards required 226.737 and 244.290 seconds, totaling 471.027 seconds (636.9
  conformers/s).
- The identical command reused both completed shards in 1.49 seconds with
  43,008 KB peak RSS and byte-identical catalog/shard manifest hashes.
- E028 produced 200/200 chemical SDFs and vdW annotations across equal 1X8B and
  8BJU lanes while preserving retrieval membership and admission order.
- Geometry QC found close receptor points in 200/200 and severe center-distance
  points in 195/200. Interpreted this as evidence that rigid retrieval overlays
  are not docking poses. Locked vdW distribution and representative chemical
  pose review before directional/torsion final-ranking integration.

## 2026-09-09 - Analyzed real E028 vdW distributions by query

- 1X8B/824 was systematically more collision-prone than 8BJU/QT9: median soft
  penalties 50.197 versus 17.337, median clashing fractions 0.582 versus 0.389,
  and median severe-atom counts 25.5 versus 15.0.
- Even the lowest-penalty representatives retained substantial clashes. The
  highest-penalty poses had 0.432 and 0.212 angstrom minimum protein/candidate
  center distances, proving physical overlap rather than marginal vdW contact.
- Locked six LOW/MEDIAN/HIGH chemical SDFs for human pose inspection. The next
  decision distinguishes terminal-group conflicts, pervasive scaffold
  penetration, and frame errors before choosing torsion-only versus full local
  docking relaxation.

## 2026-09-09 - Locked E030 after human pose review

- Human inspection confirmed all representatives occupy the correct pocket and
  ruled out a gross coordinate-frame error.
- 8BJU collisions were primarily terminal without visible macrocycle
  penetration. 1X8B LOW was terminal dominated, while MEDIAN/HIGH cores sat too
  close; HIGH appeared visually better than MEDIAN despite its larger scalar
  vdW penalty.
- Locked E030 to integrate bounded terminal torsions and 13 small rigid seeds on
  a fixed final Top-N. Baseline/min-vdW/max-directional/balanced variants remain
  separately auditable, and no collision metric alters protected admission.

## 2026-09-09 - Implemented E030 deterministic rigid micro-seeds

- Added identity, +/-0.25 angstrom Cartesian translations, and +/-5 degree
  Cartesian rotations about the candidate centroid: 13 fixed seeds total.
- Verified exact identity retention, deterministic matrices/order, centroid
  preservation under rotations, and invariant pairwise atom distances.
- Chemical geometry tests pass 5/5 and the full suite passes 120/120. This is a
  reusable primitive only; directional/vdW/torsion Top-N variant integration is
  still pending and no production rank has changed.

## 2026-09-09 - Pivoted from WEE1 pose repair to general fast E031

- The user correctly identified that the project requires a general fast 3D
  search followed by later stages, not WEE1-specific pose optimization inside
  retrieval.
- Quantified the risk: 13 rigid seeds times a bounded two-torsion beam can reach
  about 403 evaluations per pose, or over three million states for the real
  7,704-member refinement set.
- Stopped E030 before any batch caller was added; current production latency and
  rankings are unchanged. Retained the tested seed primitive only for a future
  docking adapter.
- Locked E031: one directional evaluation on each existing rigid pose, <=10%
  latency-overhead target, single-pose vdW only at capped docking handoff, and
  no torsion/rigid micro-search in the general retrieval lane.

## 2026-09-09 - Clarified the final 3D retrieval match degree

- The user specified that layered retrieval should culminate in one score for
  reproduction of co-crystal ligand key interactions, not pose optimization.
- Defined `interaction_match_score` as weighted query-anchor coverage from a
  deterministic one-to-one typed spatial/directional feature assignment under
  each already retained rigid pose, with per-anchor assignments retained.
- Kept shape and ordinary color as earlier evidence lanes rather than choosing
  arbitrary combined weights. Current anchor-weighted atom-centered Gaussian is
  an approximation, not the requested complete interaction score.
- Recorded the generality gap: direct hydrogen-bond anchors exist; salt bridge,
  aromatic, cation-pi, hydrophobic, and metal interaction extractors still need
  explicit target-independent validation.

## 2026-09-09 - Implemented E031 key-interaction matching sidecar

- Implemented deterministic maximum-weight one-to-one assignment of typed,
  spatial, and signed/axial directional matches under each existing rigid pose.
- Preserved the Gaussian result byte-for-byte and wrote scores, assignments,
  per-anchor contributions, input hashes, timing, rank correlation, and Top-K
  overlap only to a separate NPZ/JSON sidecar.
- Added a Linux runner for the accepted 8BJU/QT9 and 1X8B/824 refinement sets,
  including `/usr/bin/time -v` resource reports and top-hit explanations.
- Offline suite passed 123/123. An exploratory in-memory kernel loop completed
  23,112 assignments (7,704 candidates x three poses) in 3.72 seconds on the
  local Windows environment. Real companion mmap latency and the <=10% gate
  remain to be measured on the workstation.

## 2026-09-09 - Recorded E031 workstation continuation checkpoint

- Added `to_human/resume-checkpoint-2026-09-09-e031-implemented.md` as the
  authoritative continuation handoff after code commit `f6cd02d`.
- Consolidated accepted E029/E028 evidence, human pose observations, the E030
  speed-protection pivot, E031 implementation semantics, real input/output
  paths, exact workstation command, and the locked 7.799-second QT9 gate.
- Marked E031 real dual-query timing and effectiveness analysis as pending; no
  exploratory timing or two-query WEE1 evidence was promoted to production.

## 2026-09-11 — Expanded library preprocessing requested (E032)

- User reported 365 MOL2 files / 277.08 GiB in Linux `mc_data`, 24 CPUs,
  approximately 58 GiB available RAM and 6255.46 GiB free disk. Old two-shard
  v1 and chemistry catalogs remain available. Linux HEAD reported as fb26c18.
- Added a separate batch driver and written E032 protocol. Preserve old IDs
  and data, use transactional per-file registration, immutable input inventory,
  hash-checked reuse and recoverable partial directories. No production DB is
  modified. Use bounded FAISS batches, shared all-shard training and independent
  index shards followed by one merge rather than cumulative full snapshots.
- Offline suite passed 132 tests, including 9 new safety/identity tests. The
  local test interpreter lacks RDKit/FAISS: Linux pilot, real chemistry parsing,
  index training/merge and full-scale timing remain unvalidated.
- Added `to_human/E032_LINUX_RUN.md` with explicit file transfer, pilot and
  full-run instructions. No remote execution, Git push or full-run completion.
- E031 real dual-query timing remains pending; E032 does not run query scoring,
  docking, relaxation or target-specific model prediction.

## 2026-09-14 — Diagnosed E032 MOL2 kekulization failure (E032b)

- The Linux pilot reached new source `N5_0.mol2` and aborted because strict
  `MolFromMol2Block(..., sanitize=True)` returned `None` for multiple aromatic
  macrocycles; multiprocessing only propagated the first rejected conformer and
  was not the root cause.
- Locked E032b before implementation. Added a shared loader that keeps strict
  parsing primary and permits an aromatic-graph fallback only when all RDKit
  sanitization operations other than `SANITIZE_KEKULIZE` pass.
- Applied the loader to artifact-v1 and chemical-companion workers and added
  strict/fallback counts to both manifests. Other sanitization failures remain
  fatal.
- Offline suite passed 132 tests with one RDKit-only module skipped. Separate
  RDKit 2025.09.2 checks passed strict/fallback/invalid rejection. The supplied
  representative failing macrocycle passed five conformers through both workers
  with 60 USRCAT values, 37 heavy atoms and 25 features per conformer.
- Linux retry on actual `N5_0.mol2`, multiprocessing, recovery behavior, FAISS
  work and full-scale timing remain confirmatory work.
- Pushed the protocol and implementation through `bdcb912` to `origin/main` for
  workstation pull and retry.

## 2026-09-14 — E032b Linux pilot confirmed

- User ran the patched pilot with Python 3.11.16, RDKit 2026.03.5 and FAISS
  1.14.3. Recovery preserved the previous `.N5_0.partial` before rebuilding.
- `N5_0.mol2` completed 5,744 inserted conformers, zero duplicates, in 6.7
  seconds. Both artifact-v1 and chemical-companion manifests reported exactly
  4,485 strict plus 1,259 controlled aromatic-no-kekulize records.
- The equality of both stage totals confirms no silent record loss and no
  chemistry/artifact disagreement. E032b is confirmatory for this real shard.
- Direction: continue the frozen full E032 batch. Full-library catalogs, shared
  FAISS training/merge, integrity checks, COMPLETE marker, recall calibration
  and total timing remain pending.

## 2026-09-15 — E032 conformer-name/index conflict triage
User reports split_0119.mol2 fails after 100,000 records with a ValueError for c--L-dA-Lnme-Wnme-dL-VNMe-c conf1. This is distinct from the compact primary ID collision fixed in bab0643: same name/index has different raw content hash. Real record differences remain unknown. Added explicit existing/incoming source, zero-based record index and hash diagnostics, plus a read-only comparison script covering committed peers and same-file peers after rollback. Targeted tests: 15 passed (PYTHONPATH=src); initial collection without PYTHONPATH failed, then corrected. No automatic renumbering/skipping, source mutation, database migration or remote deployment. User confirms previous source/output paths. Standalone Linux diagnostic command: D:\agent\to_human\AIDD_split_0119_diagnostic.sh. Await real output before choosing identity-policy repair; full batch remains incomplete.

## 2026-09-15 — E032 same-file label collision confirmed; optional preservation
User diagnostic confirms split_0119 zero-based records 131568 and 131625 share the molecule/conf1 label, each with 116 atoms / 118 bonds and matching ordered topology hash, but different coordinates and raw hashes. No committed peer. Implemented explicit --preserve-index-conflicts: preserve each distinct record under a separate conformer ID and negative internal index; retain original name/hash/source position, audit conflict in warnings_json, expose per-source conflict count. Existing IDs, source hashes, artifacts and default strict policy unchanged. Exact duplicates remain deduplicated; topology failures retain whole-source rollback. Downstream source-position/ID mapping reviewed. Full local suite 136 passed, 1 RDKit-only module skipped. Real Linux resumption and stereo identity QC remain pending; no claim of chemical equivalence.

## 2026-09-15 — E033 assessment runner ready
User requests acceptance plus retrieval speed, quality and candidate-reduction measurements, with GitHub handoff. Protocol frozen in 1e94c8c before local checks. Added a read-only runner that requires E032 COMPLETE, holds existing batch lock, checks generated hashes and registry/artifact/chemical identity, handles relocated old source paths via stem/hash, streams exact descriptor truth, excludes complete query molecules, measures a fixed budget/nprobe grid and exports real retained counts and cap scenarios. No filtering percentage or activity quality invented. Local suite: 147 passed, 2 skips (RDKit module and real FAISS test unavailable in current interpreter); Bash syntax passed. Real-interface test remains included for workstation environments. Isolated FAISS dependency download attempted but did not complete; no real local FAISS benchmark claimed. Synthetic end-to-end report pipeline and corrupt-data failures covered. Ready to push; actual workstation acceptance and scientific target-specific validation remain pending.

## 2026-09-15 — User checkpoint: waiting, E033 not yet validated
Latest user clarification supersedes the earlier registration-complete wording: remaining data registration is still awaited, then workstation validation will run. E032 COMPLETE is not verified; E033 has not been executed/validated on the real library. Implementation and documentation were successfully pushed through bbd7675 (including 80acb64). Saved to_human/resume-checkpoint-2026-09-15-e033-awaiting-validation.md and updated research-state.yaml. Next action is to receive completion status or report.md/report.json after scripts/run_e033_library_acceptance.sh. This turn only records progress; no new batch, benchmark or docking run.


## 2026-09-18 — E033 workstation results received
User supplied E033 report text and partial JSON: acceptance passed, 25,813,808 conformers and 8,318,351 source-grouped molecules. Eight-query calibration panel passes at candidate_budget=10000/nprobe=128: minimum recall .99, mean .997125; search p50/p95 .00635649/.00720884 seconds, search+fetch/exact-rerank median .04594543 seconds, retained molecules 8352–9834. This is protocol-aligned engineering evidence reported by the user; complete raw report and integrity hashes have not been independently inspected here. No new local benchmark or biological/pose acceptance claimed. User confirms workstation report directory /mnt/local/hand/yuzhang/aidd/e033-library-acceptance/20260918-091604. Supersedes the September 15 awaiting-results state. Next: generate query-matched WEE1 candidates on the accepted expanded index, then Gaussian/E031 measurements with matching IDs and same-run timing baseline. Existing E031 shell runner hardcodes the old 77.99-second QT9 baseline; do not apply it to new candidate sets. Snapshot: to_human/resume-checkpoint-2026-09-18-e033-panel-passed.md; summary: data/e033-user-reported-20260918.json. No business-code changes or remote runs in this checkpoint.

## 2026-09-18 — E034 implementation and local validation
Protocol d2aaa73 was committed before implementation checks. Added expanded_wee1.py and run_e034_expanded_wee1.sh: verify E033 completion/report and current batch lineage, snapshot locked 8BJU/QT9 and 1X8B/824 queries, generate query-specific descriptors, compare nprobe128/256 at budget10000 against full-corpus exact Top1000, gate both queries at .95 before Gaussian/E031. No borrowed library-panel candidates or early one-conformer cap. Completed stages use input/output hash receipts; interrupted Gaussian chunks resume, and partial-refine timing makes the <=10% gate unavailable. E031 full-call timing uses the same-run query-specific Gaussian refine wall time; legacy fixed77.99 comparison removed. Added current RUN_STATUS and process-local open-file soft-limit handling for many mmap shards. Local suite166 passed/2 dependency skips;19 new tests include actual Gaussian and E031 with synthetic mmap payloads, tampering, recall failure gating and resume. Bash syntax and JSON/YAML checks passed. RDKit/FAISS unavailable in current Windows interpreter; isolated dependency download did not complete, so real chemical preparation/FAISS and physical25.8M E034 remain workstation-pending. User requested GitHub push and workstation handoff: to_human/E034_WORKSTATION_RUN.md. E033 supplied summary also saved in tracked to_human/e033-user-reported-20260918.json. No raw library uploaded.

2026-09-18 — E034 startup correction: workstation command with --resume failed because the selected output has no protocol.json. This occurred before E034 execution, not evidence of library or recall failure. Driver now treats --resume on a nonexistent output as a fresh run, while preserving and rejecting existing directories lacking protocol.json with actionable path guidance.21 focused E034 tests passed, including fresh-start with/without --resume, actual resume, changed-configuration rejection and existing-content preservation. Previous full suite166passed/2skips remains historical; real workstation E034 still pending.

2026-09-18 — E034 blocked at E033 evidence: user reports EVALUATION_COMPLETE.json absent from confirmed report.md directory. Exact cause (wrong path, incomplete copy or interrupted finalization) unverified. Added aggregate missing-file diagnostic with clean CLI exit; no fabricated marker or bypass. Added run_e034_with_fresh_acceptance.sh to create unique new evaluation/output directories, rerun read-only E033 acceptance/calibration and enter E034 only on successful evaluation. Existing library and reports preserved.23 E034 tests passed and recovery Bash syntax passed. User workstation evaluation remains pending.

## 2026-09-18 — E034 workstation result received
User supplied completed report.md and partial report.json from /mnt/local/hand/yuzhang/aidd/e034-expanded-wee1/recheck-20260918-163847-316118. Both locked WEE1 queries have exact-USRCAT Top1000 recall1.0 at budget10000/nprobe128 and256; select128. QT9 retains6299conformers/2212molecules after refine, refine31.242s, E0314.870s, overhead15.59% fails10% gate.824 retains6791/2961, refine69.954s, E0314.421s, overhead6.32% passes. This is protocol-aligned user-reported engineering evidence, not an independent artifact/hash audit. Exploratory ranking disagreement: anchored Gaussian versus E031 rho .1801/.1919; conformer-level Top100 intersections1/2. Do not promote E031 to ranking or infer biological/pose superiority. Next reuse these artifacts for molecule-level disagreement/pose review and numerical-equivalence-preserving E031 profiling; do not rerun registration or change the gate. Details: to_human/E034_RESULT_2026-09-18.md and .json. No new computation in this result-recording turn.

## 2026-09-18 — E035 implemented; local large synthetic equivalence passed
User requested molecule-level disagreement review, equivalent E031 acceleration and large-scale validation while keeping annotation-only policy. Protocol1674d11 committed before checks. Added opt-in batched geometry and feature-only mmap reader; original Hungarian solver and reference default unchanged. Molecule comparison chooses independent Gaussian/E031 representatives with global-ID tie order and exports both when distinct, cross-scores, anchor contributions, native-frame heavy-atom SDF, crystal copy and PyMOL review script. Full runner verifies original E034 hashes, profiles reference separately, performs3alternating full-call repeats and gates exact assignments/ranks plus abs(score/anchor error)<=1e-6. Historical E034 timing denominator is labeled and requires host/CPU-count match; actual workstation timing remains pending.
Default scale samples100000distinct real library conformers across all shards and compares both query sets x3constructed centroid-aligned orientations (600000comparisons); optional1M/6M. These are constructed stress poses, not refined/docking evidence. Hash-bound512-conformer checkpoints support resume. Local synthetic100000cases passed, max error1.1102230246251565e-16, zero assignment mismatch and exactfloat32 ranks. Nonempty kernel groups speedups2.10–3.68x, excludes I/O/serialization and does not establish QT9<=3.1242s full-call target. Full local suite188passed/2dependency skips; end-to-end synthetic mmap test covers profile/review/SDF coordinates and charges/scale/resume/tamper checks. Bash syntax passed. Handoff:to_human/E035_WORKSTATION_RUN.md; metrics:to_human/E035_LOCAL_KERNEL_100K.json. Physical100k/1M and human pose review await user workstation; no Gaussian rerun or ranking mutation.

2026-09-18 — E035 final review: scale reference now uses the original full artifact/chemical readers independently and checks every sampled feature coordinate/type/direction/ID against the optimized reader before scoring. Added induced reader-corruption rejection test. Final suite189passed/2dependency skips; Bash syntax and state/metrics parse checks passed. All original E034 inputs remain read-only; default E031 engine and admission/ranking policies unchanged.

## 2026-09-18 — E035 million result and E036 scheduling implementation
User supplied E035 completion: 1,000,000 distinct conformers / 952,402 source-grouped molecules / 6,000,000 comparisons; zero score error and assignment mismatches; all global ranks exact. QT9/824 full optimized E031 medians1.447585/1.909238s, speedups2.4996/2.2462 versus same-experiment reference, historical E034 overhead4.63%/2.73%. Evidence is user-pasted, not independently rehashed remote artifacts. Molecular Top100 shared6/9; exported39/41poses still await human review. Annotation-only remains mandatory. Summary:to_human/E035_WORKSTATION_1M_RESULT.json.
E036 protocol70beeff fixed before tests. Implemented bounded process scheduling (2x workers), spawn E035 scale workers with independent mmap readers, deterministic assembly, original hash receipts, partial resume, progress/ETA and separate cumulative compute versus elapsed wall time. E035 wrapper defaults8workers/512rows with a NEW output directory. Gaussian now bounds in-flight tasks too. New fast_3d_search entry point loads the existing whole-library FAISS index once for two calibrated WEE1 queries, retains nprobe128/budget10000, uses coarse500/refine64/16workers and optimized E031; verifies all Gaussian arrays and E031 assignments/scores/ranks against E034. No changed seed cap, admission or ranking. Exact whole-library reference scans and million stress validation are absent from this query path; startup integrity/index-load and final checks are separately labeled. Fresh candidate schedules must exactly match original E034 evidence.
CONFIRMATORY local correctness:194passed/2dependency skips; real spawn mmap serial/parallel equality, partial resume/tamper rejection, Gaussian chunk-size/worker equality, bounded failure handling and complete indexed-query pipeline with fake FAISS tested. Bash syntax and diff checks passed. Actual chemistry/FAISS workstation performance remains pending; no parallel speedup claimed. Tiny synthetic correctness timings are not performance evidence. User handoff:to_human/E036_FAST_3D_RUN.md. Prior completed run-1m-v1 remains untouched; old code-hashed runs must not be resumed under new implementation. Next measure actual Gaussian latency and inspect exported discordant poses; if scheduling is insufficient, profile worker kernels before reducing search effort.

## 2026-09-18 — E037 local progress while workstation unavailable
User requested further autonomous implementation and one-shot workstation validation. Protocol1a46816 precedes tests. Found deterministic pair_alignment_seeds generated all seeds then sliced512; added opt-in max_seeds early return after exactly the same unique prefix. Gaussian reference default retained; bounded flag is hash-bound in scaled config. All candidate budgets/Top-N/seed caps/scoring and E031 annotation-only policy unchanged. Added first-chunk worker cProfile/text diagnostic in separate profile runs, never latency-eligible.
Consolidated workstation suite crosses old/fine chunk layouts with reference/bounded seed generation, two rounds in reverse order, then one optional separate profile scenario; default wrapper enables profile. 16 timed query executions plus2profiled. Each child runs existing full-index search and compares all outputs against E034. Suite preserves per-run logs, hash-bound receipts, partial-resume exclusion, failure continuation and final nonzero status on failures; no million E035 rerun. CLI: scripts/run_e037_workstation_suite.sh. Handoff:to_human/E037_ONE_SHOT_VALIDATION.md.
CONFIRMATORY local correctness:205passed/2dependency skips. Tests cover seeded/degenerate/symmetric seed prefixes, cap boundaries, full Gaussian arrays, multiprocessing/profile artifacts, orchestration order/failure/resume and exclusion of profiled/partial-reuse timings. Local synthetic12cases x3alternating generator repeats all exact; capped-case generator speedups1.27–9.69x, uncapped~0.99–1.01x. Raw timings:to_human/E037_LOCAL_SEED_BENCHMARK.json. These are synthetic generator results, NOT full-query or workstation speedups. Real chemistry/FAISS performance and exported-pose review remain pending.

## 2026-09-19 — E038 prompt/protein/AF3 interface implemented
User requested LLM prompts for protein data, AF3 and subsequent 3D workflows, confirmed OpenAI-compatible provider and AF3 already installed on workstation. Protocol1ba97f4 preceded implementation tests. Added strict JSON intent planner (fixed actions/typed dependencies/no model-supplied argv, paths or sequences), compatible chat transport with environment-only credentials and redacted HTTP errors, active-Project ownership and existing AI audit integration, hashed plan/stage receipts and explicit compute capability. Public UniProt fetch/unique gene+taxonomy resolution validates accession/species/sequence; ambiguous matches fail rather than choose. RCSB candidates/downloads reuse existing adapters; no receptor selection is inferred. AF3 preparation takes verified sequence, optional explicit1-based construct and validated CCD ligands; local profile compiles argv with shell=False, imports model/confidence/provenance, preserves failed attempts. Only one protein chain plus optional CCD ligands is currently exposed. E031 remains annotation-only. Search prompt adapter selects only calibrated WEE1 QT9/824/both; fast_3d_search now supports one-query execution without changing scientific parameters.
CONFIRMATORY local engineering validation:222passed/2dependency skips full regression; final AF3 path-type change rechecked16targeted tests. Standalone credential-free synthetic smoke passed, including resume. Mock tests exercise provider request shape/response validation, invalid actions/dependencies, secret redaction, organism/ambiguity, ownership/tampering, AF3 prepare/run argv and output import, compute blocking/resume and query adapter. No live model API, protein API, real AF3/GPU or library performance run here. Offline fixtures are labeled, not scientific results. State/JSON/Bash syntax checked. Handoff:to_human/E038_PROMPT_PROTEIN_AF3.md; summary:to_human/E038_LOCAL_VALIDATION.json. Tomorrow:run_e038_prompt_smoke.sh, E037 suite if pending, then local LLM/AF3 config and run_e038_prompt.sh. E037 passing supports calibrated engineering use, not new-target/pose/activity acceptance.


## 2026-09-20 - E039 English content and workflow skills

Translated maintained guides, script messages, examples and historical checkpoints
into English while preserving original evidence dates and scientific measurements.
The previously untracked September 15 checkpoint is retained as an English historical
record. Added repository instructions and four discoverable skills for calibrated
3D search, verified protein preparation, installed AF3, and prompt orchestration.
The planner receives compact skill/action groups; validated envelopes record routing.
English summaries/questions are requested and Han text is rejected. A repository
check covers tracked and nonignored UTF-8 content, excluding raw ignored data and
Git history. README now reflects supplied E034/E035 results and pending E037 timings.
Validation: 224 tests passed, two dependency skips; all four skill validators passed;
E033 shell syntax and whitespace checks passed. This is implementation validation,
not a new scientific experiment. Live LLM/AF3 and end-to-end optimized workstation
measurements remain pending. E031 remains annotation-only.


## 2026-09-20 - Consolidated workstation acceptance handoff

Added to_human/WORKSTATION_FULL_ACCEPTANCE.md after inspecting the current wrappers,
planner actions, receipts and E037 report fields. Distinguishes supplied E033-E035
evidence from pending E037 timings and live E038 model/protein/AF3/search execution.
Includes one combined live compute prompt, reuse checks and noncompute boundary
cases. No new workstation or scientific execution was performed.


## 2026-09-21 - E040 container AF3 and provider profiles

Implemented Apptainer execution matching the supplied SIF, weights/databases, GPU
flag and XLA setting. Native execution remains compatible. Output bind directories
are created before launch; failures preserve attempt directories; image SHA-256 is
computed once per invocation and bound to resume. Added independent GPT/DeepSeek
profiles and explicit provider selection, preserving the legacy profile environment.
Setup copies local profiles without overwriting existing files. Keys remain separate
environment variables. Confirmatory software validation: 230 passed, two dependency
skips; shell syntax, English content and updated skill metadata passed. Mock tests
cover command mapping, missing host inputs, failed/retried/reused execution, changed
image rejection and provider credential isolation. No live API or GPU run performed.
Official current DeepSeek setup uses deepseek-flash; GPT profile uses gpt-5.6-terra.
See to_human/AF3_CONTAINER_AND_LLM_PROVIDERS.md for workstation commands and sources.


## 2026-09-21 - E040 user-reported live AF3 completion

DeepSeek Flash planned live retrieval of human WEE1 P30291 (646 residues), AF3
preparation and installed Apptainer GPU prediction. All three steps completed;
RTX 5090, seed 1, five samples. Data pipeline 420.33 s (MSA 412.50 s); inference
66.27 s. These are nested stages, not an independently measured total latency.
pTM 0.51, ranking score 0.79, fraction_disordered 0.56; no pose-quality acceptance.
Confidence JSON reports seven repeated A IDs despite single-chain arrays; the
local adapter copies that JSON unchanged, so raw summary/CIF inspection is pending.
Evidence is user-pasted report/log; see to_human/E040_WORKSTATION_AF3_RESULT.json.
E037 performance, prompt-driven search, GPT and workstation reuse remain pending.


## 2026-09-21 - Live prompt-driven WEE1 search and raw AF3 follow-up

User supplied a complete live DeepSeek execution report for wee1_both in
PROMPT-4689d49ad1834f89. E031 remains annotation-only. Inner search timings and
per-query equivalence details have not been supplied, so no speedup claim is made.
Raw AF3 summary also contains a long repeated-A list while chain_ptm and pair
matrices describe one chain. This confirms the repetition is present upstream of
our confidence importer, not that the model contains multiple chains. Token-level
labeling remains a hypothesis; preserve raw evidence without silently deduplicating.
No AF3 rerun is needed solely to inspect metadata. See
 to_human/E040_WORKSTATION_SEARCH_RESULT.json. GPT and E037 timing remain pending.


## 2026-09-21 - E041 conversational workbench

Implemented loopback browser chat, persisted sessions/messages/jobs, separate model
selection, a serialized background subprocess queue, report-derived status/results,
process-tree cancellation, explicit receipt-based resume and owned-run attachment.
Scientific reports now expose running/completed steps atomically while work proceeds.
HTTP tests cover bearer authentication, origin checks and session task isolation.
Provider router uses bounded recent dialogue and task summaries; raw local files
are not automatically uploaded. Maintained UI and responses remain English.

User supplied fresh search timings: QT9 32.859 s (Gaussian 28.641, E031 1.696),
824 64.056 s (Gaussian 62.087, E031 1.938), all reported equivalence checks passed.
Startup/index loading and final checks are separate; this is not an E037 speedup
comparison. Updated the existing workstation evidence record without claiming a
new independently executed benchmark.

User reports Schrodinger/Maestro via modules and PLANTS installed. Added a read-only
installation probe; actual docking CLI/license/grid/reference inputs remain pending.
New-target query calibration, docking, rescoring and biological validation are not
marked complete by adding a chat UI. The guide records stage-by-stage acceptance.
No live model/GPU chat test was performed here. Repository checkpoints provide
continuity because this environment exposes neither /loop nor a cron tool.

E041 validation: full regression 240 passed, two dependency skips; final chat
tests 10 passed; desktop screenshot visually inspected. Scientific computations
were not repeated. Remaining live chat and docking discovery steps are documented.

### 2026-09-21 Docking installation handoff
Operator supplied Schrodinger 2025-1 paths and a PLANTS home directory. Added a separate trusted discovery profile, home expansion, permission checks and explicit-path precedence. Three local fixture tests passed. No binaries were executed; licenses, receptor/grid inputs, live docking and chat docking integration remain pending. See to_human/DOCKING_INSTALLATION.md.
Full local regression: 243 passed, 2 optional-dependency skips; English guard and Bash syntax checks passed.

## 2026-09-21 - E042 screening evidence and explicit selection

Protocol b91d479 locked the same-pose selection and immutable handoff contract.
Implemented session-bound review/select/export actions using completed search
receipts, source hashes and separately requested export. Crystal anchor provenance
is mapped to stored feature columns; legacy reconstruction must match every query
array. Reports classify HBA/HBD, retain original residue/geometry evidence, provide
per-query and union counts, and explicitly leave other interaction classes unassessed.
Selection evaluates ALL/ANY in one conformer/objective pose before molecule
aggregation and an optional explicit cap. Gaussian ranks and E031 annotations are
not rewritten. Exports preserve stable IDs and heavy-atom rigid poses; docking,
ligand preparation and biological quality remain unvalidated.

Confirmatory local regression: 258 passed, 2 optional-dependency skips. Fifteen
E042 tests passed again after a platform-independent fixture adjustment. Tests
cover same-pose versus mixed-conformer evidence, exact legacy mapping reconstruction,
zero results, stable representatives, missing/invalid conditions, cross-query
rejection, source modification failures, session isolation, and a real child-process
review/preview/empty-export chain. Nonempty SDF export uses mock library readers.
English-content and both updated skill validators passed. Live model routing,
real-library export and manual crystal/pose review remain workstation acceptance.
PLANTS executable permission success is user-reported; no docking was run here.
Guide: to_human/E042_SCREENING_TO_DOCKING_HANDOFF.md. No AF3 or library search was
repeated locally. No continuity scheduler is exposed; repository checkpoints retain state.

### 2026-09-21 E042 workstation ID representation fix
User-reported evidence review failed at retrieval/refinement molecule identity comparison. Production retrieval stores S16 byte IDs; Gaussian stores U16 text IDs. str(bytes) introduced a false mismatch. Normalize identifier arrays with strict UTF-8 decoding on read, preserving scientific arrays and genuine mismatch checks. Updated fixtures to production byte retrieval IDs, plus genuine mismatch and invalid encoding tests. Full regression: 260 passed, 2 skips; English guard and diff checks passed. Keep completed search artifacts; restart updated chat and create a fresh evidence task rather than resuming an old code-bound review. Live retry pending.

## 2026-09-21 — E043 exhaustive screening and classified query evidence

Implemented exact chunked descriptor Top-K over all catalog conformers, a native
crystal-query feature classification report (HTML/CSV/JSON), and classified
selection/export provenance. PLIP and candidate complexes are not prerequisites.
Seven interaction hypothesis categories are represented; halogen-specific matching
is explicitly unsupported by the current feature schema. Original E031 ranking
is unchanged. Added reference-ligand pocket derivation and separate PrepWizard /
LigPrep / Glide commands with a reference redocking gate before candidate docking.

Local exploratory engineering validation: 267 tests passed, 2 dependency skips.
Full-library timing, real crystal chemistry and licensed 2025-1 execution remain
pending on the workstation. No activity or candidate pose-quality acceptance.
Protocol: experiments/E043-exhaustive-contacts-docking-protocol.md.
Guide: to_human/E043_EXHAUSTIVE_CLASSIFIED_SEARCH.md.

## 2026-09-21 — E044 uncapped condition evaluation

User corrected the membership requirement: count the whole library under feature
conditions, not the Top-10,000 descriptor candidates. E043 workstation evidence
confirmed exhaustive descriptor coverage for both WEE1 queries, with 4,214/4,042
retrieved molecules and 2,254/2,902 refined molecules (union 5,138). Classification
also ran on the workstation, but these budgeted sets cannot answer full-library
condition counts.

Implemented a separate full_count chat action: persistent workers evaluate every
catalog conformer with no descriptor Top-K or Gaussian Top-N; write score/assignment/
pose chunks; merge molecule counts on disk. Existing diagnostic levels are shown
without asking the user for another pre-compute threshold. Same-pose all/any
selection streams the completed chunks and supports separately requested SDF/ID
export. Partial runs never establish full coverage; changed/corrupt chunks fail.
The algorithm still chooses one heuristic anchored-Gaussian pose per conformer,
so no exhaustive orientation/torsion or biological recall claim is made.

Exploratory engineering validation: 272 tests passed, 2 dependency skips; the
final input-validation guard passed its 5 targeted tests. English guard: 240 files.
New tests cover last-chunk hits, duplicate molecules across chunks, partial-to-full
resume, corruption, same-pose conditions, zero/nonzero export, local chat routing
and unchanged Gaussian pose computation. Licensed docking and actual uncapped
full-library throughput/counts remain untested here. Guide:
to_human/E044_FULL_LIBRARY_CONDITIONS.md.

## 2026-09-21 — E045 necessary-condition funnel after observed E044 cost

User log: QT9 completed 18 additional 2,048-row chunks between 14:21:29 and
14:28:38, or 36,864 rows / 429 s = 85.93 rows/s. Local-rate extrapolation gives
about 83 hours for remaining QT9 work, excluding the second query and aggregation;
it is not a reliable full-run ETA. E044 removed budgets but did not implement
cheap necessary-condition rejection, so it missed the intended fast-funnel design.

Added an explicit-rule funnel using type/direction-kind compatibility, injective
feature assignment and conservative rigid pair-distance bounds before the same
Gaussian pose search. All survivors are scored without Top-K/Top-N; actual passes
are deduplicated. Conditional feature counts are labeled as such. Unscored rows
have explicit masks and cannot be exported as calculated poses. Rule changes
require a new funnel; no threshold or required contact set is silently invented.
The unfiltered per-feature diagnostic question remains distinct and potentially
expensive, especially for ANY/single-feature rules.

Compatible frozen E044 chunks can be reused through the CLI, with input/scoring
code/parameter checks and a hard failure if a bound rejects a previously passing
pose. Progress now counts completed rows rather than interpreting the largest
out-of-order ID as cumulative coverage, and reports new pose work separately.

Exploratory local validation: 279 passed, 2 dependency skips. Tests cover exact
pipeline pass IDs on a small fixture, type/matching/pair-distance rejection,
random rigid-transform retention at four thresholds, boundary cases, source-task
routing, and legacy contradiction detection. Real necessary-condition survivor
rates and end-to-end speedup are pending. Guide:
to_human/E045_NECESSARY_CONDITION_FUNNEL.md.


## 2026-09-21 - E046 after a nonselective E045 workstation pilot

User stopped the run. Displayed 55,296 conformers all survived necessary bounds;
final matches were sparse. This falsifies practical pruning efficiency on that
prefix, not correctness or final specificity. Protocol was written before tests.
Implemented seed-batched float64 Gaussian overlaps with shared color kernels and
reference recomputation of near winners. Added conservative relative-direction
bounds, affinity-sized CPU pools, optional one-owner CuPy, timing diagnostics and
a session-bound /benchmark action. The pilot scores even rejected sample rows and
checks all Gaussian/E031 arrays and actual condition membership against reference.
Local five-conformer exploratory fixture: complete arrays exact; median 4.04x CPU
scoring speedup, excluding library I/O and E031. GPU unavailable locally. No full
library acceptance, CUDA performance or improved biological discrimination claim.
Guide: to_human/E046_HARDWARE_FUNNEL.md. Larger cross-conformer GPU batching and
quantified real-library pruning remain open if this pilot is insufficient.

E046 regression: 300 passed, 3 dependency/GPU skips. English guard: 250 maintained files passed. No physical CUDA device was exercised in this environment.


## 2026-09-21 - E047: change the work, not only the device

User E046 GPU report: 256 samples, exact arrays, no positive matches, every sample
survived invariant bounds. NumPy/22 median 1.586 s; CuPy/1 median 7.143 s, second
repeat 4.868 s. Real CPU speedup did not establish acceptable library scaling.
Raw supplied timings are transcribed in to_human/E046_WORKSTATION_GPU_RESULT.json.

Implemented a cheap rule-specific test over the unchanged existing seed set,
before whole-shape Gaussian evaluation. Independent anchor maxima upper-bound
any assignment; ALL remains within one seed. Skip a conformer only if all seeds
are impossible, and keep every original seed competing on survivors. Reuse seed
construction/decoded records. Preserve full ID coverage and explicit unscored masks.
No invented contacts, tighter threshold, Top-K or altered winner policy.

The benchmark now supplements spread IDs with saved positive/near-threshold rows,
reports current reference positives, checks skipped negatives and exact surviving
arrays, and separates seed, feasibility and Gaussian timing. GPU remains optional.
This does not yet remove linear seed-generation work. A generic index of the same
nonselective invariant bounds was not built: it would not justify a speed claim.
Real-library benefit of E047 remains pending a workstation pilot.

E047 local regression: 313 passed, 3 dependency/GPU skips; English guard 255 files; skill validation passed. Passing actual E031 fixtures, same-pose ALL, threshold boundaries, preserved original seed competition, full ID coverage with skipped Gaussian rows, and false-rejection detection are covered. Workstation speedup and whole-library totals remain unmeasured.


## 2026-09-21 - E048 joint coarse eligibility

User clarified that shape, chemistry and crystal-derived anchors must jointly
constrain membership before seed generation. E047 user-reported 288-row pilot
retained all invariant rows; 199 reached Gaussian after seed feasibility.
NumPy22 medians were 1.526214 s enabled and 1.580669 s disabled; current
reference positives were zero and scalability_gate was false. This is not
acceptance of the pre-seed funnel.

Implemented explicit atom-ratio, principal-extent and typed-feature coverage
criteria before existing anchor necessary bounds. These add eligibility rules;
they do not preserve the old anchor-only pass set by definition. Threaded policy
through chat preview, full funnel and reference comparison. Added /coarse for
10000 spread IDs plus saved evidence, with no seed/Gaussian computation.
Pocket-derived features can be requested; receptor collision/occupancy checks
are not implemented. All reported new selectivity remains workstation-pending.

E048 local validation: 324 passed, 3 dependency/GPU skips; English guard passed
259 maintained files; repository search skill validated. Tests include joint
predicate rigid invariance, explicit policy propagation, pre-seed rejection,
full ID accounting, pose-free audit and session-bound /coarse routing.


## 2026-09-21 - E049 survivor validation

User-reported E048 coarse audit retained 288/10032 in 3.674 s without seeds.
The smaller hardware panel retained 41/288 before seeds, then zero for Gaussian;
reference positives zero, scalability gate false. Added a terminal entry point
that reconstructs the larger panel, saves IDs, compares every coarse survivor
and sampled rejects, and reports per-anchor failures and synthetic query controls.
This is not whole-library acceptance or an independent chemical positive control.

E049 local regression: 328 passed, 3 skipped. English guard passed 263 files.
Synthetic identity/Gaussian checks and detection of lost real-panel positives
are fixture-tested; workstation control results remain pending.


## 2026-09-21 - user-authorized overnight full-library validation

Added an orchestration entry point for pre-run controls, uncapped current-rule
whole-library evaluation, complete molecule enumeration, and sampled final-hit
reference reproduction. Uses 22 CPU workers, bounded chunks, no Top-K/Top-N.
Current scope is the explicitly selected QT9 query; no 824 rule is invented.
Full-library stage resumes verified chunks with unchanged code/inputs/settings.
Wall-time completion overnight remains unmeasured.

Overnight orchestration regression: 331 passed, 3 skipped; English guard 266 files.
Tests reject incomplete coverage, keep uncapped execution parameters, report
missing positives, and deterministically sample across final result lists.


## 2026-09-21 - workstation run started; overnight handoff

User confirms the full-library validation is now running and appears much faster.
This is user-reported progress, not observed completion or a measured speedup.
Expected code is 1533883; actual PID, exact output directory, command and current
coverage have not been supplied. Do not launch a duplicate or update running code.
Tomorrow inspect suite-report.json, full-library/progress.json, run.log and final
coverage/counts. Keep scientific acceptance pending until reports are available.
Handoff: to_human/20260921_OVERNIGHT_HANDOFF.md.

## 2026-09-22 - correct preselection funnel validation scope

User supplied completed E049 ALL run: 25,813,808 conformers, 8,318,351 molecules,
647,638 coarse survivors, zero Gaussian-evaluated conformers and final hits.
Subsequent user-run controls passed crystal direction consistency and three
synthetic rigid recoveries; twelve selected conformers had no four-anchor hits
under original or expanded tested seeds. These are user-reported artifacts,
not locally ingested complete workstation data.

User clarified the intended funnel: preserve one-or-more anchor matches, record
which anchors and pose-specific combinations, cluster, then allow human selection
before downstream docking/affinity work. Early four-anchor ALL was incorrect for
this scope. User specifically requires substantial computational reduction before
human selection, with measured retention and speed. Stop ALL-zero-hit diagnostic
expansion. E050 protocol records a corrected pilot and scale-up plan; production
code and workstation execution are unchanged. Corrected funnel acceptance pending.

## 2026-09-22 - direct full-library E050 entry delivered

User explicitly requested direct full-library execution and rejected small-data
pilot prerequisites. Added independent preselection_full entry and launcher;
existing E049 scoring modules are unchanged. New protocol uses ANY over the
selected anchor IDs with the recorded coarse rules and 0.5 threshold, processes
every catalog conformer, retains a real pose for every observed anchor bitmask,
merges molecules without synthesizing simultaneous contacts, and groups retained
members by exact non-stereochemical Murcko scaffold. Acyclic connectivity and
explicit singleton failures avoid silent member loss. This grouping is not 3D
similarity clustering. No docking/affinity stage is launched.

Output includes sealed resumable chunks, all-stage conformer/molecule counts,
worker compute and wall-time scopes, full candidate pose/member tables and a
SQLite index. Family/output locks prevent duplicate same-family jobs. Actual
workstation selectivity, wall time and candidate quality remain unmeasured.
Regression: 338 passed, 4 skipped; the new RDKit scaffold reconstruction test
was skipped because RDKit is unavailable locally. Fixture orchestration covers
full coverage, resume, ANY routing and pose evidence retention. English guard
passed 272 maintained files. Run guide: to_human/E050_FULL_PRESELECTION.md.
