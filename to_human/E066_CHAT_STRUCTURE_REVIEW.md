# E066: Chat structure review, SMILES complexes and desktop PyMOL

Use the [workstation checklist and remaining-work register](E066_WORKSTATION_CHECKLIST.md)
for the ordered 2026-09-25 tests, pass criteria and acceptance record.

Implementation and local adapter tests are available. Real AF3 inference, live
provider routing and PyMOL GUI acceptance must be tested on the workstation.
This update does not run a new full-library search or modify old scientific outputs.
Local validation: 515 tests passed, 4 skipped; English guard and JavaScript syntax
checks passed. The update is maintained on feature/structure-guided-chat. Prefer
`git pull --ff-only origin feature/structure-guided-chat` on that branch after
checking that the working tree is clean. Workstation deployment remains separate.

For the previously supplied offline patch only: copy `E066_CHAT_UPDATE.patch` to a separate checkout at base 34457db (or a compatible
checkout). From its repository root run `git apply --check E066_CHAT_UPDATE.patch`,
then `git apply E066_CHAT_UPDATE.patch`. Stop if the check fails; do not force it or
apply into an active scientific job's checkout. Do not apply the patch after pulling
the same changes. The patch includes the previously
prepared offline E065 replay utility; it does not activate production clustering.

## Workstation setup

Keep an active scientific job on its original checkout. Test this version in a
separate checkout, or update after that job finishes. Restart the Chat server to
load new code; retain its existing storage/project configuration to see old tasks.
Do not resume old sealed compute receipts under changed code. Reviewing completed
results is read-only; new predictions/searches use new task/output directories.

Use the existing AIDD environment and AF3 runtime profile. Add this key to the
existing runtime JSON, retaining its AF3, library and other settings:

```json
{"pymol": {"executable": "/absolute/path/to/pymol"}}
```

Find the installed desktop executable with `command -v pymol`. PyMOL runs on the
Chat server host, outside the AF3 container. Start Chat from a working desktop
session with the appropriate DISPLAY, or use your existing remote desktop.
Opening Chat in a browser on another computer does not launch that computer's
PyMOL. A headless SSH session needs an explicitly configured graphical session.
The bridge reports startup failures and writes `viewers/SESSION/pymol.log` in
Chat storage. No network PyMOL RPC server or arbitrary model-generated Python is used.

Launch with the same configuration as before:

```bash
conda activate aidd-workstation
export AIDD_RUNTIME_PROFILE=/absolute/path/to/your-existing-runtime.json
export AIDD_LLM_CONFIG_DIR=/absolute/path/to/your-existing-llm-config
bash scripts/run_chat_agent.sh --allow-compute
```

RDKit must be installed in the AIDD environment for strict SMILES validation.
The existing AF3 installation, databases, weights and native/container profile
are reused. Do not reinstall AF3 for this update.

## Experimental structures and diverse same-pocket ligands

Example English prompts (replace task IDs and reference identifiers as needed):

```text
Retrieve human MDM2 and find ligand-bound experimental structures.
Compare the ligands and propose structurally diverse references.
Use 3JZK, chain A, ligand YIN at residue 1 as the reference pocket. Build the same-pocket consensus from the diversity task TASK_ID.
Open the diverse reference complexes from task TASK_ID in PyMOL.
Open all admitted complexes from task TASK_ID in PyMOL.
```

Use the existing diversity/consensus workflow; its actual admission report
determines which complexes share the pocket. PyMOL applies the saved consensus
alignment to displayed structures. A diversity-only task loads original coordinates;
explicit display alignment is available but does not establish pocket equivalence.
The viewer supports up to 100 structures at once. The `all_admitted` option applies
to a consensus; a diversity report contains only its proposed references.

For deterministic opening, enter `/pymol TASK_ID`. The response lists object IDs
such as `v001` and their experimental query IDs. A queued operation is not yet
successful: check `/pymol_status` for a completed or failed receipt.

## AF3 with your own ligand SMILES

Use a verified protein accession and an explicit construct when appropriate.
Replace the bracketed SMILES placeholder with the actual complete string, including
stereochemistry and formal charges. Do not send the placeholder literally.

```text
Retrieve human MDM2 Q00987. Prepare an AF3 complex using residues 25-109 and one ligand with this exact SMILES: [PASTE_SMILES_HERE]. Use seed 1. Prepare inputs only.
Retrieve human MDM2 Q00987. Run AF3 for residues 25-109 with one ligand using this exact SMILES: [PASTE_SMILES_HERE]. Use seed 1 and the configured AF3 installation.
```

The construct above is an example for a defined domain, not automatic selection
of a biologically optimal construct. Input preparation stores `af3-input.json`,
`ligand-chemistry.json` and protein/input provenance in the task directory. Review
the plan and input before expensive inference when using a new ligand or construct.

The adapter supports one protein chain and up to eight total CCD/SMILES entities.
SMILES is copied unchanged after strict parsing; canonical SMILES is diagnostic
only. Disconnected fragments must be separate entities. Unspecified tetrahedral
stereocenters are reported, not filled in. No silent protonation/tautomer enumeration,
custom covalent bonds, multiple protein chains or affinity prediction is added.
Successful RDKit parsing does not guarantee AF3 conformer generation succeeds.
An inference failure remains failed, with logs and prepared inputs retained.

## Read and interpret prediction results

```text
Show the results of task TASK_ID.
Explain the confidence metrics of AF3 task TASK_ID.
Open the predicted complex from task TASK_ID in PyMOL.
```

Deterministic commands: `/results TASK_ID`, `/confidence TASK_ID`, `/pymol TASK_ID`.
The report exposes the actual model and `confidence-analysis.json` paths. Analysis
contains only metrics present in the output; missing values are not fabricated.

| Metric | Interpretation |
| --- | --- |
| ipTM | Confidence in relative placement across chains; not binding affinity or proof of binding. |
| pTM | Confidence in the overall predicted structure, not specifically the ligand interface. |
| chain_pair_iptm | Off-diagonal interface confidence; diagonal within-chain pTM. Use recorded chain IDs to interpret matrix order. |
| pLDDT | Local per-atom confidence, 0-100, stored in AF3 model B-factor fields. Not an experimental temperature factor. |
| PAE | Predicted aligned error in angstroms; lower means more confident relative placement. |
| ranking_score | AF3 sample ranking, combining confidence, disorder and clash penalty; not a probability. |
| has_clash | AF3 clash flag; inspect the geometry and ligand chemistry. |

The summary adapter explains pLDDT/PAE but does not currently derive full per-atom
or per-token plots. Those arrays remain in the AF3 output files. No universal
ligand-confidence cutoff or biological acceptance is asserted.

## Control the open PyMOL session

Natural-language examples:

```text
Show the protein as a cartoon and the ligands as sticks.
Show polar contacts between each ligand and its own protein within 3.6 angstroms.
Show the pocket surface with opacity 0.4.
Color the ligand in v001 magenta and zoom to its pocket.
Label protein residue 96 in chain A of v001.
Hide v002 and rotate the view 45 degrees around y.
Save a PNG snapshot and then save the PyMOL session.
```

If a provider misroutes a compound instruction, use one operation per message.
These direct commands bypass provider interpretation:

```text
/view {"operation":"cartoon","target":"protein"}
/view {"operation":"sticks","target":"ligand"}
/view {"operation":"polar_contacts","cutoff":3.6}
/view {"operation":"surface","target":"pocket"}
/view {"operation":"transparency","target":"pocket","opacity":0.4}
/view {"operation":"color","objects":["v001"],"target":"ligand","color":"magenta"}
/view {"operation":"zoom","objects":["v001"],"target":"pocket"}
/view {"operation":"label_residues","objects":["v001"],"target":"protein","chain":"A","residue":"96"}
/view {"operation":"hide","objects":["v002"]}
/view {"operation":"rotate","axis":"y","angle":45}
/view {"operation":"snapshot"}
/view {"operation":"save_session"}
/pymol_status
```

Other operations: `contacts` (heavy-atom proximity, default 4.5 angstroms),
`hide_contacts`, `hide_surface`, `hide_labels`, `show`, `background`, and
`align` with an explicit reference object ID. `show`/`hide` enable/disable entire
objects; rotation changes the view. `align` performs display-only protein CA
alignment. The pocket display selection is protein residues within 5 angstroms
of the selected ligand. Opening a new catalogue clears objects managed by this
Chat bridge, but does not delete unrelated user-created PyMOL objects.

Each complex is handled separately: no inter-complex contact lines are generated.
Polar lines are PyMOL geometry/typing candidates, not verified hydrogen bonds.
Protonation and chemistry are not independently validated. Aromatic, hydrophobic,
halogen, metal and water-mediated interactions are not all inferred by this command.
Contact CSV/JSON includes atom names, chains, residues and distances. The Chat
sidebar shows the latest receipt and download buttons for its PNG, PSE, CSV or
JSON artifacts. Earlier artifact paths remain in their operation receipts.
After restarting PyMOL, reopen the task before issuing display controls.

## Inspect existing crystal contact evidence

```text
List the recorded interactions in consensus task TASK_ID.
For consensus task TASK_ID, show contacts involving residue 96, including atom names, distances and the recorded interaction classes.
Show up to 200 recorded contact pairs for query 3JZK:YIN:A:1 from task TASK_ID.
```

Use `/interactions TASK_ID` for the first 50 pairs, class counts, policy notes,
available query IDs and full ledger path. Natural-language filtering supports
query ID, residue, contact class and up to 500 rows. Author and canonical residue
numbers are retained; inspect both when mapping to the displayed experimental
structure. The returned class definitions distinguish proximity from confirmed
bonds. These are recorded consensus contacts, not newly validated interactions.
For AF3/PDB-only tasks without a consensus ledger, request PyMOL contacts instead.

## Tomorrow's acceptance sequence

1. Open one completed consensus task; compare displayed ligand alignment with
   its recorded reference. Check selected object labels and ligand chain/residue.
2. Display polar contacts; inspect several distances manually and download CSV.
   Change color/surface, hide an object, export PNG and reopen the saved PSE.
3. Prepare a known workstation-tested ligand SMILES and verified construct.
   Compare the generated input with the earlier successful standalone AF3 input.
4. Run that small known complex; compare reported confidence with AF3's actual
   summary JSON; open its final model in PyMOL and inspect the ligand.
5. Run the separate [seed equivalence and width audit](E064_SEED_BUDGET_AUDIT.md)
   on completed search data. Do not infer production speed from unit tests.
6. Only after molecular equivalence passes, request a fresh Chat budget search
   with explicit `seed_search` settings; do not change a resumed run's protocol.

Example opt-in search prompt, using an existing completed consensus design task:

```text
From design task TASK_ID, start a fresh budget search with 20 workers, chunk_conformers 4048, retrieval_molecules 1000000, export_molecules 100000, and seed_search {"backend":"batched","max_pair_seeds":512}. Export original MOL2 records and molecule names.
```

Check the generated plan's parameter names and sealed protocol. Default seed
backend and cap remain unchanged; this update adds explicit Chat access, not
automatic adoption. Current retrieval is the existing USRCAT/FAISS route, not a
new validated pharmacophore index. Target-specific recall remains unmeasured.

## References

- [Official AF3 input specification](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)
- [Official AF3 output specification](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)
- [PyMOL Python querying API source](https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/querying.py)
