# Target Structure and Ligand Similarity Workflow

## Structure sources

All experimental and predicted structures enter one candidate registry. Source
types include PDB experimental, AlphaFold Database, local AlphaFold 3,
AlphaFold Server import, AlphaFold 2, and user-supplied structures. Raw files
are immutable Artifacts; prepared receptors are separate derived Artifacts.

## Multiple-structure review

A Campaign may freeze a named comparison set with at least two candidates, one
alignment reference, chain assignments, pocket residues, bound-ligand IDs, and
a scientific rationale. Selection of a comparison set does not select the final
receptor.

```text
target resolution
-> RCSB candidate search and metadata
-> raw mmCIF acquisition
-> comparison-set selection
-> pocket alignment and deviation analysis
-> PyMOL review
-> final receptor choice
-> receptor preparation
```

Pocket analysis aligns matched C-alpha coordinates by Kabsch superposition and
reports global pocket RMSD, per-residue displacement, coverage, and missing
residues. Ligands are displayed and compared separately because ligand identity
and pose can differ even when protein pockets are similar.

RCSB nonpolymer entities are stored in two views. `nonpolymer_ids` preserves
every reported component, while `ligand_ids` is a conservative candidate-ligand
view that excludes recognized water, solvents, ions, buffers, and
crystallization additives. Every exclusion and its reason remains in
`excluded_nonpolymer_components`; unknown IDs stay in the ligand view for human
review.

## PyMOL checkpoint

The Agent writes a `.pml` file into the active Project. The client launches
PyMOL locally; a headless desktop or cluster node is not required. The script
loads candidates, aligns them to the reference, displays the selected pocket,
shows organic ligands, groups objects, and orients the view. Final selection
requires a recorded user rationale.

## Fast ligand similarity

Similarity search is an indexed library service, not an on-demand scan of MOL2
text files.

### 2D

1. Standardize molecule identity and canonical SMILES.
2. Build Morgan fingerprints once per library/version.
3. Persist fingerprint blocks with molecule IDs and index parameters.
4. Use vectorized/Bulk Tanimoto search.
5. Return ranked unique molecules with provenance.

### 3D

1. Build USRCAT descriptors for every registered conformer.
2. Use the descriptor matrix for a fast top-N prefilter.
3. Rerank only the shortlist with aligned shape/color similarity.
4. Aggregate conformer hits to unique molecules while retaining the best
   conformer and score.

The index key includes the library hash, molecule standardization version,
fingerprint or descriptor parameters, and RDKit version. A changed library or
parameter set creates a new immutable index.

RDKit supplies chemistry operations; Gemmi/BioPython supply structure parsing;
PyMOL supplies interactive visualization. They are optional deployment
environments rather than mandatory control-plane imports.
## Predicted structures in a Campaign

AlphaFold 3 outputs are discovery candidates, not automatic replacements for
experimental receptors. Import a completed output directory with its construct
and model version:

```bash
aidd-agent import-alphafold3-result --db /path/to/aidd.sqlite3 \
  --user USR-... --campaign CAM-... --output-dir /path/to/af3/output \
  --input /path/to/input.json --construct target_kinase_299_569 \
  --model-version source-commit-or-release --chain A \
  --manifest-output /path/to/prediction-manifest.json
```

`campaign-status` then reports `predicted_candidates` beside the RCSB
`candidates`. Import does not change the scientific selection. After comparison
and human review, either candidate type can be frozen by internal candidate ID:

```bash
aidd-agent select-receptor --db /path/to/aidd.sqlite3 --user USR-... \
  --campaign CAM-... --candidate PRD-... --rationale "Reviewed rationale"
```

The lock snapshots the exact candidate metadata and structure SHA-256.

## Review the complete receptor ensemble

Use one PDB only as the alignment reference; all Campaign PDB and predicted
candidates remain peer members:

```bash
aidd-agent create-receptor-ensemble --db /path/to/aidd.sqlite3 \
  --user USR-... --campaign CAM-... --name all-target-structures \
  --reference-pdb 8BJU --pocket-residues 320,337,463 \
  --exclude-pdb "9TG7=WEE1 is only a 12-residue degron peptide, not the kinase domain" \
  --rationale "Observe receptor and ligand-pocket diversity before selection"
```

Generate a PyMOL script after all experimental mmCIF files are present under
the Project `inputs/structures` directory:

```bash
aidd-agent prepare-receptor-ensemble-pymol --db /path/to/aidd.sqlite3 \
  --user USR-... --ensemble ENS-... --project-root /path/to/project \
  --output /path/to/project/target/reviews/ENS-.../review.pml
```

The script aligns proteins with `cealign`, displays only classified ligand IDs,
and selects the 6 Å protein neighborhood of each ligand. Predicted apo models do
not receive invented ligands.

## Ask an LLM to review receptor eligibility

Write JSON describing the intended receptor purpose and required domain. This
creates an audited, provider-neutral prompt/evidence packet; it does not call a
model and does not exclude a structure.

```bash
aidd-agent create-receptor-eligibility-review \
  --db /path/to/aidd.sqlite3 --user USR-... --campaign CAM-... \
  --requirements-json receptor-requirements.json \
  --computed-evidence-json receptor-computed-evidence.json
```

A later provider adapter must return exactly one `include`, `exclude`, or
`manual_review` decision for every candidate. Import it with
`import-ai-recommendation`. Unknown evidence, missing candidates, duplicates,
and `requires_human_review=false` are rejected. Accepted output is advisory; a
human explicitly creates final ensemble exclusions and the receptor lock.
