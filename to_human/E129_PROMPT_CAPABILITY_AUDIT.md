# Prompt capability audit: full-library search and AF3 dimers

Date: 2026-11-11.

## Supported paths

| Request | Chat route | Execution |
| --- | --- | --- |
| Protein monomer, optional ligand | workflow -> af3_prepare -> af3_run | Configured native/container AF3 |
| Homodimer, optional ligand | One sequence with protein_copies=2, or two identical entries | Two protein entities A/B in one complex |
| Heterodimer, optional ligand | Two literal sequences or two verified protein dependencies | Two distinct protein entities A/B in one complex |
| Arbitrary-target full-library search | target_search prepare/start/run | Target query -> guided_funnel -> preselection_full |
| Search using a dimer receptor | Explicit receptor_chains A/B plus posed reference | Both selected chains contribute to receptor geometry |

Ligands are optional. With two proteins and one ligand, the ligand is entity C.
Each ligand entry belongs to the same complex; independent compound predictions
require separate prepare/run pairs. Input chain IDs must be distinguished from
potentially malformed output summary labels (see E128).

## Prompts

Replace placeholders with complete sequences or exact absolute paths on the
Agent host. Avoid claiming a seed count alone uniquely identifies a seed list.

```text
Run AlphaFold 3 for a homodimer: exactly two identical protein chains in ONE complex.
Use this exact sequence for both chains, without trimming or mutations.
Protein sequence:
<FULL_SEQUENCE>

Do not include any ligand.
Model seeds: [47, 48, 49, 50, 51]
Prepare the input JSON, run the configured AF3 runtime, and keep all outputs locally.
Report the task ID, input JSON, model CIF and confidence paths.
```

```text
Run AlphaFold 3 for a heterodimer in ONE complex: one copy of sequence A
and one copy of sequence B. Preserve both complete sequences without mutations.
Sequence A: <FULL_SEQUENCE_A>
Sequence B: <FULL_SEQUENCE_B>
Do not include any ligand.
Model seeds: [47, 48, 49, 50, 51]
Run the prediction and report the task ID and saved input/model/confidence paths.
```

For one ligand in either complex, replace the no-ligand sentence with:

```text
Include exactly one ligand in this complex with this exact SMILES: <SMILES>
Preserve its chemistry and formal charges.
```

```text
Search the entire configured local molecular library.
Use receptor /absolute/path/receptor.cif and receptor chains A and B.
Use posed reference ligand /absolute/path/reference.sdf in the receptor coordinate frame.
Name this query DIMER_SITE_SEARCH. Use target_search to prepare and start the search.
Report the applied screening policy, task ID, full_coverage, processed conformer
count and saved candidate/output paths. Keep candidate data on the local host.
Do not run docking.
```

Select only chain A for a monomer receptor. The reference must define the actual
site in receptor coordinates; a SMILES or centered SDF alone does not define a
binding pocket. A predicted AF3 complex can provide a structural hypothesis, but
automatically extracting its ligand and starting a search is not established by
this audit. Provide the prepared receptor/reference files explicitly.

## Configuration and limits

Restart the Agent after updating, retaining storage/runtime/provider settings.
Execution needs `--allow-compute`; AF3 needs a working `af3_profile`, while search
needs `search.batch`, `workers` and `refine_chunk`. External molecular input files
must be under the Project or configured `--input-root`. GPU memory/runtime can
limit large complexes; two full CIP2A copies double the protein residue count.

Full coverage means all catalog conformers are considered under the recorded
coarse filters and finite pose search. It does not mean all possible conformations
or poses are searched, all binders retained, or predicted affinity established.
The engine sets full_coverage only after processed count equals catalog count.

Generic full-library candidate storage is local. Automatic generic MOL2 quota
handoff remains pending; budget-workflow exports are separate. The optimized AF3
shared-MSA batch runner remains CLI-only. Do not claim either as completed Chat
functionality. A request containing only a target name still needs a selected
structure and site/reference before this generic search can execute.

## Audit evidence and fix

Four real Chat-worker subprocess checks exercised homodimer/heterodimer with and
without ligand, retaining exact sequences, IDs and seeds. Both the planner and AF3
scientific runner were synthetic. Existing tests covered query preparation,
same-pose contacts, compute gates, container commands and full-catalog scans.

The audit found that target_search continued invoking the language model after
successful queue submission. It now returns a deterministic task receipt and
guards against another target_search mutation while an operation is pending in
the same agent turn. This does not deduplicate deliberately repeated future turns.

No live provider, GPU dimer inference or production-library timing/recall run was
performed here. The earlier user-returned CIP2A single-protein/ligand report is
evidence for that execution only, not live dimer acceptance.

Validation: 111 distinct targeted tests passed across two runs; English-content
guard passed 739 maintained files; diff check passed.
