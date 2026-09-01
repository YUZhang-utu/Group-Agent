# Findings

## Current understanding

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
