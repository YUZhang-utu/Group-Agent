# E103: EquiScore evaluation of the completed PLANTS panel

Protocol before execution, 2026-10-05. User selected EquiScore and requested an
installation script because it is not installed on the workstation. No inference
is yet claimed. Use an isolated Linux environment, leaving AIDD dependencies intact.

Freeze official repository commit 8b2a9289cf7d181fa49de6ac6712260e8c500c4a and
save_model_screen.pt SHA256 d4367bb73686b2363e238abb778fab55e2924458ec1ced561072bd82f711695d.
The official Screening.py uses softmax class 1 and descending ordering. Preserve
this direction and score meaning; it is not calibrated experimental affinity.

Run a technical pilot on 32 deterministically selected conformer IDs (hash order,
seed 20261005), each against all available receptors. Selection is independent of
ChemPLP score. Validate source seals, same-frame receptor heavy atoms, source/pose
ligand graph identity, original coordinates, checkpoint compatibility and explicit
per-pair missingness. Any failed pilot pair blocks the full run until corrected.
This gate establishes executable compatibility only, not scientific model accuracy.

Full evaluation reuses all saved conformer/receptor poses, including candidates
outside ChemPLP's Top-10. Deduplicate computation by CID/receptor. Score one frozen
pose per pair; no redocking or minimization. After inference, independently choose
the highest-scoring in-block conformer per molecule and the ten highest-scoring
molecules per block/receptor. Retain Top-5 views and compare ranks with ChemPLP.
Do not average model units. Retain failures; incomplete panels remain unranked.

Official preprocessing uses whole residues within 8 Angstrom of each ligand and
RDKit/ProLIF features. Our runner follows that feature path while rejecting errors
instead of silently removing metals or dropping samples. Model weights load with
strict key matching. Graph exceptions and inference failures remain explicit.

Primary outputs: per-pair EquiScore and original pose identity, per-block Top-10
means/ranks, ChemPLP-versus-EquiScore rank comparisons and selected-candidate overlap.
Without measured activity labels, these remain exploratory score comparisons, not
EF1, binding affinity validation or proof of improved partition quality.

Engineering checks precede workstation use: higher-is-better ranking, molecule
deduplication, source/hash binding, unknown/duplicate/missing result rejection,
pilot/full distinction, interrupted execution receipts and no docking dispatch.
