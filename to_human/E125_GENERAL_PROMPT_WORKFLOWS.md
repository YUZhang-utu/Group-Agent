# General prompt-driven AF3 and full-library screening

Date: 2026-10-29. This extends E124 from report visibility to generic execution.
The scientific search does not dispatch by KRAS/WEE1 target name. Block-routing
development remains deferred. All UI and maintained workflow text is English.

## Update and start

On the workstation, in the existing Group-Agent checkout and active environment:

```bash
git pull --ff-only origin feature/structure-guided-chat
python -m pip install -e .
```

Restart Chat using the existing storage root, runtime and provider configuration.
Keep those settings to preserve conversations and the configured AF3/library paths.
For explicitly supplied local receptor/reference files, append an allowed input root:

```bash
--input-root "$HOME/aidd-inputs"
```

An illustrative full command follows. Replace the runtime path with the existing
local configuration; do not replace a working storage directory with a new one:

```bash
python -m aidd_agent.chat_web \
  --storage-root "$HOME/aidd-chat" \
  --runtime /absolute/path/to/runtime.local.json \
  --input-root "$HOME/aidd-inputs" \
  --allow-compute
```

The runtime uses the existing `af3_profile` plus a target-independent `search`
section containing `batch`, `workers`, and `refine_chunk`. No new installation,
database copy or API key in Git is required. AF3 must run on the host with the
installed native/container runtime and GPU; file paths refer to the Chat host.
Opening the browser on another computer does not move execution to that browser.

## AF3 through natural-language prompts

Replace angle-bracket placeholders with actual sequences/SMILES:

```text
Predict this protein as a monomer using AF3 and five seeds starting at 1:
<SEQUENCE_A>
Run the prediction and show its saved confidence and model path.
```

```text
Predict a homodimer of this exact sequence with seeds 7, 19 and 31:
<SEQUENCE_A>
Include one ligand with SMILES <SMILES>. Run AF3.
```

```text
Prepare an AF3 heteromer input: two copies of sequence A and one copy of sequence B.
A: <SEQUENCE_A>
B: <SEQUENCE_B>
Use ten seeds starting at 1. Add one GDP CCD ligand. Prepare only; do not run inference.
```

The planner copies supplied sequences literally, normalizing only case/whitespace.
It must not invent or silently fetch a different sequence. Supplied sequences are
recorded as user-provided, rather than claiming UniProt validation. Alternatively,
existing UniProt lookup steps may provide one or several verified sequences.

Typed AF3 fields:

| Field | Meaning |
| --- | --- |
| `protein_sequences` | List of supplied sequences, one entry per protein chain; repeat entries for stoichiometry |
| `protein_step` | Existing single verified-protein dependency |
| `protein_steps` | List of verified-protein dependencies; allows heteromers/repeated chains |
| `protein_copies` | 1..8 identical copies when only one sequence source is supplied |
| `seeds` | 1..100 explicit distinct integers |
| `seed_count`, `seed_start` | Alternative deterministic seed range; default start is 1 |
| `ligand_smiles`, `ligand_ccd` | Optional ligand entities in this one complex |

Choose exactly one protein-source mode and one seed mode. Protein and ligand
chains have distinct IDs, at most 26 total; at most eight ligand entities.
For multiple proteins, provide the intended constructs; a single global start/end
interval cannot ambiguously truncate different chains. Unsupported amino-acid
symbols, disconnected SMILES and custom covalent bonds require explicit preparation.

Multiple SMILES in ONE input mean co-ligands in ONE complex. For independent ligand
predictions, explicitly ask for separate inputs/prepare-run pairs. The optimized
large-batch runner still runs through CLI and can reuse shared protein features;
Chat separate prepare/run pairs do not currently share that MSA cache automatically.
Input validation does not guarantee GPU memory capacity for an arbitrary complex.

## Arbitrary receptor and site through prompts

```text
Use receptor /home/y/yuzhang/aidd-inputs/enzyme.pdb, protein chain A,
and posed reference ligand /home/y/yuzhang/aidd-inputs/reference.sdf.
Call this query ENZYME_SITE1. Prepare the search query and show supported contacts.
Do not scan the library yet.
```

```text
Run a full-library search for the prepared query task <TASK_ID>.
Require contact with chain A residue 123 and chain A residue 167 in the same pose.
Show the actual policy, progress and output files.
```

An explicit receptor selection is used as the user's structural hypothesis; it
does not falsely claim a completed receptor-cohort review. If receptor selection
is still needed, the existing structure/pocket-state survey workflow remains available.
For a PDB ID, specify the receptor chain and reference ligand's CCD/author chain/
author residue identity. The adapter retrieves that PDB and CCD; ambiguous ligand
instances are rejected instead of choosing the first one. There is no target-name
whitelist. Predicted receptor/reference complexes may be supplied as files, but
their binding pose is a hypothesis, not a measured co-crystal.

The reference must already lie in the receptor coordinate frame. An unbound,
centered SDF does not define a binding site and is not automatically docked/aligned.
PDB/mmCIF receptors and single-molecule SDF/MOL2 references are supported. Multiple
receptor models require explicit selection before import. Protein atom alternate
locations use highest occupancy, then alphabetical altloc ID for ties.

## Deterministic Chat command

This uses the same executor without asking the LLM to translate the request:

```text
/target_search {"operation":"start","name":"ENZYME_SITE1","receptor":"/home/y/yuzhang/aidd-inputs/enzyme.pdb","reference_ligand":"/home/y/yuzhang/aidd-inputs/reference.sdf","receptor_chains":["A"],"required_contacts":[{"chain":"A","residue":"123"},{"chain":"A","residue":"167"}]}
```

Use `operation: prepare` to stop after query generation. Later:

```text
/target_search {"operation":"run","task_id":"RETURNED_QUERY_TASK_ID","required_contacts":[{"chain":"A","residue":"123","residue_name":"HIS"}]}
```

For a public structure, substitute `pdb_id` for `receptor` and `ligand_instance`
for `reference_ligand`, for example this already used KRAS reference (one example,
not an implementation restriction):

```text
/target_search {"operation":"prepare","name":"KRAS_REFERENCE","pdb_id":"7RPZ","receptor_chains":["A"],"ligand_instance":{"ccd":"6IC","chain":"A","residue":"202"},"keep_residues":[{"chain":"A","residue":"201","name":"GDP"},{"chain":"A","residue":"203","name":"MG"}]}
```

`keep_residues` retains specifically identified cofactors/ions in receptor exclusion.
Other hetero residues and water are excluded. Missing requested cofactors fail
preparation. Nonstandard protein residues are retained when recognized as amino
acids, but their side-chain contact templates may be incomplete.

## What a full-library run means

The adapter seals imported inputs, builds ligand chemical features and a real
receptor coordinate model, checks the reference against the requested policy, and
calls the existing guided full-catalog engine. It considers every catalog conformer
without a candidate Top-K/Top-N membership budget. Existing chunk receipts and
code/configuration checks control resume. Status/results are existing Chat task
operations; they must not submit another scan.

For each requested residue, compatible observed reference features form an ANY
group; all requested residue groups must pass in the same pose. Actual candidate
heavy-atom proximity to every requested residue is also checked in that pose
(default 4.5 A, configurable per contact as `distance`). `feature_class` can narrow
a group to a reported HBD/HBA/hydrophobic/aromatic interaction class. `residue_name`
can require the expected amino-acid identity. This is not hydrogen-bond angle/energy
validation and does not simulate internal torsions or receptor flexibility.

Without `required_contacts`, the broad exploratory default requires ANY observed
reference feature. It does not silently require every contacted residue. More
than 20 selected feature anchors requires narrowing the rule explicitly.

Default generic policy is recorded in the design report:

```json
{
  "minimum_score": 0.5,
  "coarse_constraints": {
    "heavy_atom_ratio": [0.1, 10.0],
    "maximum_extent_distance": 10.0,
    "minimum_feature_coverage": 0.01
  },
  "pocket": {
    "minimum_heavy_atom_distance": 1.2,
    "maximum_clashing_fraction": 0.0
  }
}
```

These broad coarse bounds avoid importing the old restrictive whole-ligand
coverage threshold as a universal default. They are not target-calibrated quality
criteria. There is no hard whole-ligand Gaussian shape/color cutoff here. Gaussian
scores rank surviving poses; the distance floor is a gross-overlap check, not
validated steric energy or cavity occupancy. A failed reference control blocks
the run rather than automatically relaxing the policy. Override the recorded
policy explicitly when warranted; do not promise a chosen number of survivors.

Output includes the prepared query/design and the existing full-search report,
candidate database, pose and molecule/scaffold summaries. This adapter does not
automatically run docking or export a fixed million-molecule MOL2 quota. Coordinate
export/docking is a separate downstream decision using the resulting candidates.

## Validation boundary

Offline tests exercise literal sequences, heteromer composition, seed expansion,
AF3 command/result dispatch, generic receptor inputs, cofactor identity, native
reference controls, same-pose contact enforcement, file seals and compute gating.
Existing Chat, full-library, guided-engine and container tests are included.
Synthetic service/dispatch fixtures are explicitly not live GPU or full-library
benchmarks. Actual throughput, target recall, model inference and docking quality
still require school-runtime acceptance and target-specific controls.

Release check: 155 offline tests passed in the isolated delivery checkout. The
English-content guard, JavaScript syntax, diff checks and Chat CLI entrypoint also
passed. Live-provider natural-language reliability is not established by fixtures.
