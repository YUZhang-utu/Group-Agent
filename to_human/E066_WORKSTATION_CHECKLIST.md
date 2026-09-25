# Workstation acceptance checklist and remaining work

Recorded: 2026-09-24. Intended test date: 2026-09-25.
Status: NOT RUN on the workstation. Local baseline: 515 passed, 4 skipped.
This is a test plan, not evidence that the desktop integration already works.

## 0. Prepare the test environment

- Update the feature/structure-guided-chat branch following the
  [setup guide](E066_CHAT_STRUCTURE_REVIEW.md). Use the previously supplied offline
  patch only as an alternative; do not apply both. An old checkout without this
  update will not have the new Chat commands.
- Use a separate checkout if a scientific job is still running. Keep old sealed
  outputs; do not resume old compute receipts using changed code.
- In the same graphical desktop session, run `command -v pymol`, then launch that
  executable manually once. Record its path and version. Close this manual window;
  the Chat bridge launches its own window and does not attach to an arbitrary
  already-running PyMOL instance.
- Add `pymol.executable` to the existing runtime JSON, preserving AF3 settings.
  Restart Chat using the existing storage and runtime/provider configuration.
- Use one completed same-pocket consensus task, called CONSENSUS_TASK below.
  If absent from the conversation, attach its existing plan through the UI's
  existing-run import field. Do not rerun the full library to test visualization.
- Prepare the exact protein construct and ligand SMILES that succeeded in your
  earlier standalone workstation AF3 test. These are the inference test fixtures.

PyMOL appears on the Chat server's desktop. It is not an embedded browser viewer
and does not launch on a different computer merely because that computer opens Chat.

## 1. Prove the direct bridge works before testing natural language

Enter `/pymol CONSENSUS_TASK`, then `/pymol_status` after the operation finishes.
Pass: a PyMOL window appears; connected is true; the operation receipt is complete;
the catalogue lists actual query labels and object IDs. Queued alone is not a pass.
For failures, save the response and `viewers/SESSION/pymol.log` under Chat storage.

Check the correct ligand chain/residue and saved alignment for at least two
complexes. Compare against the reference in the consensus report. A diversity-only
task is not a substitute here: it initially loads structures in their original frame.

Run each command separately and wait for a complete receipt:

| Command | Pass condition |
| --- | --- |
| `/view {"operation":"cartoon","target":"protein"}` | Protein cartoon is visible. |
| `/view {"operation":"sticks","target":"ligand"}` | Intended ligand atoms appear as sticks. |
| `/view {"operation":"color","objects":["v001"],"target":"ligand","color":"magenta"}` | Only the selected ligand changes color. |
| `/view {"operation":"zoom","objects":["v001"],"target":"pocket"}` | Camera focuses on that pocket. |
| `/view {"operation":"surface","objects":["v001"],"target":"pocket"}` | Pocket surface appears. |
| `/view {"operation":"transparency","objects":["v001"],"target":"pocket","opacity":0.4}` | Surface becomes translucent. |
| `/view {"operation":"rotate","axis":"y","angle":45}` | View rotates; source files remain unchanged. |
| `/view {"operation":"hide","objects":["v002"]}` | v002 disappears; skip if only one object exists. |
| `/view {"operation":"show","objects":["v002"]}` | v002 reappears. |

## 2. Verify contacts and downloaded files

```text
/view {"operation":"polar_contacts","objects":["v001"],"cutoff":3.6}
/pymol_status
```

Pass: the receipt completes with a contact count and JSON artifact. CSV is produced
when pairs exist. Inspect at least three returned atom pairs where available:
object, chain, residue, atom name and displayed distance must match the table.
No contact may connect atoms from different complex objects. Zero pairs is a valid
geometric result, not automatically a failure; use a known-contact fixture to test
line rendering. Polar candidates are not chemically validated hydrogen bonds.

Download CSV/JSON immediately from the latest-operation card. Then run:

```text
/view {"operation":"snapshot"}
/pymol_status
```

Download the PNG and verify it opens and matches the displayed view. Then run:

```text
/view {"operation":"save_session"}
/pymol_status
```

Download the PSE and reopen it manually in PyMOL; check objects and representations.
The sidebar shows the latest operation's files, so download before the next command
or retain the earlier receipt's artifact path.

## 3. Test live English routing

Use the actual object IDs returned by the catalogue. Send one action per message:

```text
Color the ligand in v001 yellow.
Zoom to the pocket of v001.
Show polar contacts for v001 within 3.6 angstroms.
Save a PNG snapshot.
```

Pass: each natural-language request performs the same intended action as its direct
command, with a complete receipt and visible effect. Record the provider/model.
If direct commands pass but English prompts fail, classify it as a routing failure,
not a PyMOL failure. Do not mark the entire Chat feature accepted using direct
commands alone. Test other configured providers separately if you will demo them.

## 4. Inspect recorded crystal evidence

```text
/interactions CONSENSUS_TASK
For consensus task CONSENSUS_TASK, show contacts involving residue 96, including atom names, distances and recorded contact classes.
```

Pass: query IDs, counts and selected pairs match the task's contacts.json ledger;
an absent residue returns no matches rather than invented interactions. Check author
versus canonical numbering before comparing with PyMOL labels. A full-ledger path
and truncation flag must be available when the table exceeds the display limit.
Use a residue known to exist in your fixture if 96 is absent.

## 5. Prepare, then run the known SMILES complex

Replace every placeholder with the previously tested values:

```text
Retrieve protein ACCESSION. Prepare an AF3 complex using residues START-END and one ligand with this exact SMILES: SMILES_STRING. Use seed 1. Prepare inputs only.
```

Pass: af3-input.json contains the correct verified construct and unchanged SMILES,
including chirality and formal charges. ligand-chemistry.json records the parse
and stereo diagnostics. This preparation-only task must not launch inference.
Compare the input against your earlier successful standalone input, not just its name.

Then submit a new task with the same explicit values:

```text
Retrieve protein ACCESSION. Run AF3 using residues START-END and one ligand with this exact SMILES: SMILES_STRING. Use seed 1 and the configured installation.
```

Call the returned ID AF3_TASK. Pass only when real inference completes and the
final model and summary exist. Prepared/queued/blocked are not inference success.
Save the execution log on failure; do not replace the installation merely because
the first Chat run fails. Inference runtime is not included in a fixed test deadline.

## 6. Interpret and display the AF3 result

```text
/results AF3_TASK
/confidence AF3_TASK
Explain the confidence metrics of AF3 task AF3_TASK.
/pymol AF3_TASK
/pymol_status
```

Pass: reported ipTM/pTM/ranking/clash values match the actual summary file. Only
available metrics are displayed; explanations do not claim affinity or activity.
The displayed structure must be the final model from that task, with the expected
protein and ligand. Repeat ligand coloring, pocket zoom, polar contacts and PNG
export. A low-confidence prediction can still pass software integration testing.

## 7. Validate the acceleration separately

Follow [E064](E064_SEED_BUDGET_AUDIT.md) using the completed search directory.
First require reference versus batched-512 exact payload/ranking equivalence, then
inspect measured timings. Compare generated caps and survivor targets for ranking
stability; do not adopt smaller budgets from speed alone. If equivalence fails,
retain the reference backend and save the failing sample/report.

After this audit, a fresh small Chat budget task can verify that explicitly supplied
workers, chunk_conformers and seed_search settings appear in the plan/protocol.
Do not start another full-library run merely to prove that Chat forwards parameters.

## Remaining work: distinguish validation from missing implementation

| Item | Current status | Next step |
| --- | --- | --- |
| Workstation deployment | Not performed here | Apply compatible patch/configuration and restart Chat. |
| Desktop PyMOL bridge/display/export | Implemented; local adapter tests passed | Execute steps 1-2 and record real GUI results. |
| Live English prompt routing | Implemented; live provider acceptance pending | Execute step 3 and retain failed prompts. |
| SMILES AF3 execution | Input adapter implemented; new real inference not run | Execute steps 5-6 using known successful chemistry. |
| Crystal contact inspection | Recorded proximity ledger and filters implemented | Check against actual ledger in step 4. |
| General arbitrary PyMOL commands | Not exposed | Extend the structured operation set when a concrete need arises. |
| Browser-embedded interactive 3D | Not implemented | Current visualization is the independent desktop window. |
| Automatic connection to any pre-existing PyMOL window | Not implemented | Use the bridge-launched window. |
| Complete validated interaction profiling | Not implemented | Add chemistry/geometry validation for required interaction classes; proximity is not a substitute. |
| Per-atom pLDDT plots, PAE heatmaps, automatic confidence coloring | Not implemented in this Chat update | Add explicit output parsing and visualization if required. |
| Multi-protein AF3, custom CCD/covalent bonds, DNA/RNA | Not exposed by current Chat adapter | Separate typed input extensions. |
| Automatic ligand protonation/tautomer/stereo enumeration | Not implemented | Current adapter preserves supplied chemistry. |
| Scientific pose/affinity/activity validation | Not established | Requires independent structural/biological controls. |
| Batched seed speedup on the workstation | Local implementation tested; real gain unmeasured | Run E064 audit. |
| New pharmacophore ANN route and target-specific recall | Not connected/validated by this update | Current budget retrieval remains USRCAT/FAISS. |
| Repair of all 139 historical Kekulize failures | Not delivered in this update | Source audit and chemistry-preserving repair remain separate; five examples do not diagnose all failures. |
| Macrocycle capacity clustering and incremental assignment | Design only | Implement after Chat acceptance; see macrocycle design. |
| Live top-tail adaptive block scheduling | Not implemented | Existing utilities are offline replay/statistics only. |

## Acceptance record to fill tomorrow

Do not change NOT RUN to PASS without observed results. Save this filled record
with task IDs, reports and screenshots; append a summary to research_log.md.

```text
Date / workstation:
Code base and patch hash:
Python / RDKit / PyMOL / AF3 versions:
Provider / model:
Chat storage / runtime profile (no secrets):
Consensus task ID / plan path:
AF3 task ID / plan path:
Direct PyMOL control: NOT RUN
Contacts and manual distance checks: NOT RUN
PNG and reopened PSE: NOT RUN
English routing: NOT RUN
Crystal ledger lookup: NOT RUN
Literal SMILES input comparison: NOT RUN
Real AF3 execution: NOT RUN
Confidence values and predicted-model display: NOT RUN
Seed audit equivalence / timing: NOT RUN
Failed prompt / exact error / log path:
Artifacts retained:
Decision / next fix:
```
