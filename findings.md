# Findings

## Current understanding

- Automatic ligand preparation must preserve crystal instance identity while
  sourcing connectivity separately from authoritative CCD definitions; mmCIF
  coordinates alone are not a safe bond-order source.
- The production MOL2 build can scan each large source file once in registry
  order, avoiding retention of all 299,999 raw MOL2 blocks in memory.
- The complete local macrocycle library can be sanitized and converted to real
  USRCAT on this laptop in about 6.2 minutes single-process with zero failures;
  parsing is not a blocker at the current 299,999-conformer scale.
- For 8BJU/QT9, broad FAISS recall is controlled mainly by IVF list coverage,
  not m20-vs-m30 PQ resolution. m20 remains viable for L1 pending more queries.
- Random access to selected chemistry, not Gaussian alignment compute, dominated
  the first L2 prototype: 7.3 seconds combined scoring versus 85.6 seconds wall
  time due to a full MOL2 rescan. Offset-addressed local binaries are justified.
- RDKit BaseFeatures matching, not MOL2 parsing, USRCAT, or PMI, dominates the
  richer artifact build (24.61 of 25.53 profiled seconds per 1,000 conformers).
  Feature topology is reusable across conformers: the first 100 three-conformer
  molecules had identical family/atom membership in all conformers, while their
  feature coordinates still require per-conformer recomputation.
- Molecule-level template reuse plus local multiprocessing makes rich 3D
  preprocessing practical at the current scale: 24 workers sustained 737.5
  conformers/s on a 30,000-conformer real sample with zero failures. The next
  constraint to measure is ordered IPC/binary writing, not feature computation.
- The full parallel writer sustained about 759 conformers/s including local
  binary I/O and checksums, completing all 299,999 conformers in 6.59 minutes.
  Int16 coordinates/features preserved source positions within the expected
  0.00501-angstrom half-bin error and USRCAT remained bit-identical.
- Offset-addressed mmap makes candidate materialization negligible at this
  scale: all coordinates/features for 10,000 candidates loaded in about 91 ms
  when sorted by file offset, versus the prior 85-second MOL2-rescan workflow.
- The library mapping is genuinely one-time and appendable: immutable artifact
  shards can be added to a trained FAISS IVF-PQ index with stable global int64
  IDs. Adding each 150K-conformer shard took about 0.3 seconds after chemistry
  preprocessing; final QT9 retrieval retained 100% of exact top-1,000 and
  99.69% of exact top-10,000 in a 100K candidate pool at about 14 ms.
- Partial three-dimensional motif admission can also reuse the immutable local
  artifacts: type-pair/distance-bin postings are built independently per shard
  from stored feature coordinates. A new co-crystal compiles to a small hashed
  invariant query and never triggers library MOL2 or feature recomputation.
- Loose, balanced, and strict pharmacophore evidence are nested views of one
  candidate channel. Their union with FAISS/L1 is recall-additive: anchors may
  admit additional local matches but cannot remove a baseline candidate.
- The exact Gaussian L2 reference now has an explicit mathematical boundary:
  identical cross/self approximation order, typed color, query-as-operand-A
  Tversky direction, parent-level alternative-point weight conservation, and
  candidate-to-query objective-specific transforms. Synthetic validation
  establishes kernel correctness but says nothing yet about real enrichment.
- Stable artifact IDs now feed that reference kernel directly through read-only
  shard maps. Candidate order and membership are invariant, while shape-only,
  atom-centered unweighted, and anchor-weighted objectives retain independent
  optimal poses and full overlap primitives. Catalog relocation is tolerated by
  resolving a missing recorded path against the catalog's local shard name.
- Artifact schema v1 cannot support honest projected-color or terminal-torsion
  claims: it lacks directional candidate projections and molecular topology.
  Those require a versioned companion artifact rather than inference from
  feature centers.
- Real E022 timing identified the Gaussian bottleneck as single-core repeated
  pose computation rather than storage: pair mode used 219.70 seconds/1,000,
  PCA-only used 16.59 seconds/1,000, both at about one logical CPU and roughly
  55 MB RSS. Therefore expensive pair alignment belongs after all-candidate PCA
  scoring, on the union of independent objective Top-N sets.
- The staged implementation now caches rigid-invariant self-overlaps, distributes
  contiguous ID chunks across processes, and commits hash-validated chunks
  atomically. Complete coarse scores remain evidence even though only Top-N
  unions receive expensive refinement; selection is compute allocation rather
  than retroactive candidate admission.
- Real staging changed the operational conclusion: all 299,999 conformers can
  receive PCA Gaussian scores in about 41 seconds and the 7,704-member Top-5,000
  objective union can receive pair refinement in about 78 seconds. The full
  Gaussian path is therefore approximately two minutes on the workstation,
  not the hours implied by the original serial reference.
- Refine CPU utilization was 743% because 7,704 candidates at chunk size 1,000
  create only eight chunks; only eight workers can be occupied. Future runs can
  use approximately 500-candidate chunks when maximizing 16-worker refinement
  utilization matters, although the present runtime is already operationally
  small.
- Anchor assignment captures observed interactions, not proven necessities.
  Anchors, projected sites, and feature counts are therefore reversible ranking
  evidence only. Broad anchor-independent recall defines admission, while
  atom-centered/projected/unweighted/anchored scores and per-anchor coverage
  remain stored for later reranking.
- A single transform is insufficient when atom-centered and projected color
  participate in optimization: each named objective must retain its own pose,
  and all reported scores must state the pose under which they were evaluated.
  Raw per-anchor cross/self overlaps plus a query-level anchor-definition
  snapshot make later normalization and human-readable reranking reproducible.
- For QT9, all three distance-supported polar contacts are also almost fully
  buried by atom-level SASA, so two hydrogen-independent observations agree.
  Reduce is absent locally; angles and Asn/Gln/His orientation remain unresolved
  instead of being inferred from missing hydrogens.

- The initial macrocycle library may contain many MOL2 molecule records per
  file; file count is not molecule count.
- Conformer identity must be modeled separately from molecule identity and later
  docking-pose identity.
- Stable internal IDs and content hashes are prerequisites for resumable
  screening and avoiding repeated computation.
- A representative existing multi-record MOL2 file parsed successfully with the
  standard-library streaming implementation.
- The complete `macrocycles-v1` library is registered locally under
  `LIB-AFA68EE6888C`: 100,000 unique molecules and 299,999 conformers, with no
  invalid input records.
- Naming-based conformer grouping produced the expected three conformers for
  99,999 molecules. The single exception lacks conf0 in the source and is safe
  to retain with conf1/conf2 unless a replacement conformer is later supplied.

## Constraints

- Physical laboratory automation requires a stricter authority boundary than
  computational orchestration. An LLM may propose an experiment but cannot
  issue free-form robot or instrument commands. Hardware execution requires a
  typed/versioned protocol, deterministic compilation, capability and safety
  validation, simulation or dry-run, and explicit human release.
- `completed` equipment execution is not a scientific result. Raw-data
  registration, assay QC, deterministic analysis, and adjudication are separate
  states, and invalid experiments must not be converted into inactive labels.
- Closed-loop learning requires exact molecule, batch, sample, container/well,
  protocol, instrument, measurement, and model lineage. Plate position or a
  filename alone is never sufficient identity.

- Source structures and experimental data remain local.
- The LLM controls workflow decisions but does not become the scientific compute
  engine or data store.
- Names are provenance/display fields, not permanent database identifiers.
- Campaign state now gates target binding, PDB candidate review, and final
  structure freezing.
- Resolution alone is not a sufficient receptor-selection rule. Real selection
  must also consider construct/mutations, binding-site completeness, relevant
  ligand/state, alternate conformations, and missing residues or atoms.
- The correct continuity boundary is the task, not the software repository.
  Shared compute artifacts can be reused by hash, but scientific decisions,
  parameters, outputs, approvals, and reports require explicit task ownership.
- An active Project is a hard precondition for Campaign and PDB operations;
  identifying only a Campaign ID is not sufficient user context.
- The operational context is `user -> active Project -> active Task ->
  Campaign`. New Campaigns must store both Project and Task identifiers, and
  Campaign/PDB operations reject mismatched active context.
- Optional methods should be inferred from actual Runs, not represented by a
  rigid pipeline or empty placeholder rows. Scientifically meaningful skips are
  decisions, not Runs.
- A uniform Run/Artifact manifest is sufficient to preserve provenance across
  local models, docking, rescoring, Slurm jobs, and MD while keeping large files
  out of SQLite.
- Daily cluster authentication is compatible with automation when credentials
  remain outside the Agent: the desktop exports a deterministic Bundle, an
  authenticated client submits it, and a content-bound receipt updates the Run.
- OpenSSH is the common transport for macOS, Linux, and modern Windows; PuTTY is
  a Windows compatibility adapter. Neither transport is allowed to disable
  host-key verification or embed passwords in Project metadata.
- A target must support a comparison set, not only one selected PDB ID. Protein
  pocket geometry, residue coverage, ligand identity, and ligand pose are
  separate comparison dimensions before final receptor selection.
- Fast ligand search requires immutable precomputed indices. Morgan/Tanimoto is
  the primary 2D route; USRCAT prefilter plus exact shape reranking avoids an
  expensive all-against-all 3D alignment.
- AlphaFold must be treated as an external compute backend, not a Python package
  inside the general agent environment. Profiles reference an independently
  installed runner, databases, and model parameters, and record explicit local
  acknowledgement of the applicable parameter terms.
- Source-control and scientific-storage boundaries are now enforceable: Git
  carries small reproducible code/configuration artifacts, while the desktop or
  cluster retains structures, libraries, indices, predictions, and credentials.
- The first real university deployment uses WEE1 (`P30291`) and has registered
  25 RCSB X-ray candidates at a maximum resolution of 3.0 angstrom. Visual
  inspection must remain distinct from the irreversible Campaign selection
  checkpoint.
- Real RCSB GraphQL responses may encode optional lists such as
  `nonpolymer_entities` as JSON null. The adapter now normalizes these values to
  empty lists instead of aborting the complete candidate search.
- The university workstation can download a selected Campaign candidate and
  open its mmCIF structure in PyMOL from the terminal. This visual review does
  not change the Campaign's formal receptor selection state.
- The current validated interface is still deterministic CLI orchestration.
  Codex has not yet been connected as the natural-language controller that
  resolves active context and invokes these tools from prompts.
- Campaign structure review is now executable end to end from public PDB IDs:
  a user can freeze a multi-PDB comparison definition, inspect it, and generate
  a Project-scoped PyMOL script without exposing registry candidate IDs or
  changing the Campaign's final receptor-selection state.
- RCSB `nonpolymer_entities` is not synonymous with bound inhibitors. Candidate
  metadata now distinguishes all nonpolymers from likely binding ligands and
  preserves excluded solvents, ions, buffers, and crystallization additives
  with reasons. Unknown components remain candidates for conservative human
  review.
- The durable AI boundary is recommendation-first rather than direct mutation:
  prompts and evidence produce provider-identified structured recommendations,
  while deterministic tools execute accepted actions. Experimental raw data
  and measurements remain outside model context; explicitly classified public
  metadata and computed summaries can support AI reasoning.
- Target prediction should not encode a permanent single "best model." The
  adapter treats AlphaFold DB as a public monomer baseline and AlphaFold 3,
  Boltz-2, and Chai-1 as versioned prediction backends. Selection must depend
  on the preparation purpose, licensing, hardware, and prospective QC.
- AlphaFold 3 is now validated as a real external discovery backend on the
  university RTX 5090 workstation. Ubiquitin and the WEE1 kinase-domain example
  completed end to end; WEE1 inference itself took about 25 seconds, while the
  roughly 47-minute MSA stage exposed NFS database I/O as the dominant cost.
- A prediction is not a selected receptor. Prediction outputs now enter a
  Campaign as immutable, hashed candidates with model, construct, confidence,
  input, and runtime provenance. The same human selection checkpoint can lock
  either an experimental PDB candidate or a predicted candidate.
- Receptor review should be ensemble-first. A named experimental structure such
  as 8BJU is an alignment reference, not an implicit winner; all experimental
  candidates and predictions remain visible, with filtered ligand identities
  and local ligand-pocket selections compared before recommendation.
- UniProt membership is not receptor eligibility. RCSB 9TG7 contains WEE1 only
  as a 12-residue degron peptide bound to beta-TrCP, so it must remain searchable
  provenance while being excluded—with rationale—from a WEE1 kinase ensemble.
- LLM receptor screening must be evidence-bounded and advisory. The interface
  requires exactly one cited decision per Campaign candidate and mandatory
  human review, so fragments, off-domain complexes, and wrong-chain structures
  can be flagged without silently modifying the Campaign.
- A PDB CCD identifier is not a unique query ligand. Campaign ligand identity
  now includes the parent structure, chain, residue, altloc, standardized
  chemistry, and optional crystal conformer. Query choice is a separate,
  immutable human lock; an LLM may recommend but cannot create it.
- The registered MOL2 library proves molecule/conformer provenance but does not
  by itself provide standardized SMILES. Real Morgan/USRCAT production indices
  therefore require an explicit chemistry preprocessing artifact rather than
  inferred or fabricated representations.

## Scalable retrieval execution

- Broad retrieval and precise reranking are complementary stored products, not
  mutually exclusive modes. Broad shard evidence remains auditable while only
  a fixed schedule receives Gaussian work.
- Billion-scale feasibility requires query-time work to depend primarily on
  shard count and configured Top-K budgets. It does not justify scoring every
  loose pharmacophore hit or writing one monolithic billion-row result.
- Baseline recall is protected structurally: FAISS IDs occupy the first schedule
  prefix and pharmacophore/anchor channels only add candidates and provenance.
- Slim coarse Gaussian chunks preserve the information needed for exact Top-N
  selection at 20 logical bytes per conformer; expensive pose detail belongs in
  the much smaller refinement product.
- Artifact v1 is sufficient for fast coordinate/feature retrieval but is not a
  substitute for the original chemistry source: it cannot authoritatively
  reconstruct atomic identity, bond order, stereochemistry, projected feature
  direction, or terminal torsions. E029 therefore consumes each original MOL2
  shard once and persists a relocatable global-ID-keyed companion.
- Physical MOL2 relocation must not mutate registered provenance. The E029
  `--source SHARD_NAME=/path` override changes only where bytes are read; every
  source record is still checked against the registry and artifact identity.
- E029 shard partials are atomic-build diagnostics, not resumable checkpoints.
  The first workstation retry exposed that source preflight previously created
  an empty partial before rejecting a missing relocated path. Source existence,
  whole-file SHA-256, and registry/artifact identity are now checked before
  partial creation; an old partial must be preserved outside the output root
  and the affected shard rebuilt.
- The Windows SQLite registry used during the original E019 build is not a
  production workstation dependency. Exact whole-file source SHA-256 proves
  the relocated Linux MOL2 bytes are identical, while artifact v1 already
  locks ordered global/conformer/molecule IDs and row shapes. E029 now uses
  those immutable inputs directly; a matching registry is optional redundant
  validation, and an unrelated active workstation registry must not be used.
- Existing 299,999-conformer throughput supports only an extrapolated
  100-million-molecule build estimate. The architecture is fixed-budget at
  query time, but 10M/100M physical scale, recall, cache, and recovery gates
  remain mandatory before production latency is reported.

## Multi-cocrystal pre-docking synthesis

- Co-crystal queries represent complementary receptor conformations and ligand
  chemotypes. Their correct aggregate is a site-scoped union, not an
  intersection and not an averaged query geometry.
- Conformer collapse must remain objective-specific: shape, unweighted color
  and anchored color can legitimately select different poses of one molecule.
- Cross-query raw Gaussian scores are not commensurate. Within-query molecule
  ranks and RRF preserve ordering under monotonic score rescaling.
- Query-protected admission lanes retain conformation-specific chemical space;
  consensus support is an additional reversible priority signal.
- A unique admitted molecule can require several docking tasks because every
  supporting co-crystal receptor is a distinct experimental hypothesis.
- Real 8BJU/QT9 and 1X8B/824 validation confirms substantial complementarity:
  only 955 of 12,834 union molecules are common to both queries (7.44%), while
  the unique sets are balanced at 5,991 and 5,888. Multi-cocrystal retrieval
  therefore needs a protected union; intersection-only selection would discard
  more than 92% of the observed candidate space.
- Of 5,000 admitted molecules, 834 support both receptor hypotheses and account
  for the extra tasks in the 5,834-task queue. Shared support is a priority
  signal, but unique-query molecules remain scientifically necessary because
  the two deposited ligand/receptor conformations explore distinct regions.
- Real E027 confirms that complementarity persists at early ranks: the two
  queries share zero Top-100 molecules, six Top-500 and 28 Top-1,000. The rank
  correlation across all 955 shared molecules is -0.0371, so shared-query
  support must not be interpreted as a stronger common raw score or affinity.
- Pre-docking geometry review can reuse artifact v1 without touching source
  MOL2: applying the retained candidate-to-query transform yields a heavy-atom
  point cloud in the crystal frame, where query coverage and conservative
  protein-distance collisions can be measured. These annotations do not alter
  retrieval membership or admission order.
- Artifact v1 cannot produce an honest SDF because it omits atomic identity,
  bonds and topology. E028 therefore labels its generic-element PDBs as point
  clouds only. A once-built v2 companion is the required boundary before
  chemistry-aware export, projected color or terminal-torsion refinement.
- E029 preserves the one-time-library boundary with a separate stable-ID
  companion rather than changing accepted v1 shards. One sequential MOL2 pass
  materializes elements, charges, aromatic/chiral flags, bonds, feature atom
  membership, directional vectors and bounded terminal-torsion move sets.
- Directional comparison must distinguish signed polar vectors from axial
  aromatic normals: reversing a donor/acceptor direction removes agreement,
  whereas reversing a ring normal represents the same plane. Translation is
  never applied to directions.
- Element-aware vdW penetration and limited torsion refinement are fixed-budget
  post-refinement operations. Zero angle remains an eligible beam state, so
  local flexibility cannot silently replace a better rigid baseline. These
  components remain reversible evidence until real redocking/enrichment
  calibrates any combined score or cutoff.
