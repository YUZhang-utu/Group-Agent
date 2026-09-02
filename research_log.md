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
