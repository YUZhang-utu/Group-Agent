# General-target pocket consensus and full-library chat test

This branch implements the E053 nonpolymer, single-admitted-receptor-state path.
MDM2 is the first real structural test, not a hard-coded target. The old E050 run
has been stopped by the user; its files are not inputs to this new workflow.

## Update and start

Use the existing scientific environment and existing runtime / LLM profiles.
From a clean checkout of the feature branch:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git status --short
git fetch origin feature/structure-guided-chat
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
python -c 'import rdkit, gemmi, Bio, numpy, scipy, requests; print("Dependencies available")'
bash scripts/run_chat_agent.sh \
  --storage-root /mnt/local/hand/yuzhang/aidd/e053-chat-workspace \
  --runtime "${AIDD_RUNTIME_PROFILE:?Set the existing runtime profile path}" \
  --llm-config-dir "${AIDD_LLM_CONFIG_DIR:?Set the existing LLM profile directory}" \
  --port 8766 --allow-compute
```

If the feature branch is already checked out in another worktree, run the pull
and start commands there. Do not reset or discard local edits. For an isolated
checkout, `git worktree add -b e053-test ../Group-Agent-e053 origin/feature/structure-guided-chat`
creates a new local testing branch. Subsequent updates there use the same
`git pull --ff-only origin feature/structure-guided-chat` command.

Open `http://localhost:8766`. Startup alone does not launch screening. All
scientific steps below are chat requests; there is no Python or JSON editing.

## Chat sequence

1. Request: "For human MDM2, accession Q00987, census experimental PDB structures
   by resolution and compare bound ligand diversity. Use 5C5A as the reference
   PDB. Do not run library screening yet."
2. After completion, request: "Build the pocket consensus using reference ligand
   5C5A:NUT:A:201 and target protein author chain A. Verify the target, organism,
   assembly, residue correspondence and aligned pocket. Keep shape templates
   independent of the full admitted interaction cohort; propose eight diverse
   templates."
3. Inspect the report's reference, admission reasons and `review.html`. An
   ambiguous reference produces explicit choices rather than silently choosing
   the first chain. Then send `/recommend`.
4. Review the cited recommendation. Explain edits in natural language, for
   example: "Make this contact optional with weight 0.4, keep either of these
   two contacts, and retain these three shape templates." Use actual displayed
   IDs or unambiguous descriptions. The adopted design must show the resolved
   IDs. Changing templates does not replace the consensus anchors.
5. Send `/adopt` to accept the recommendation, or submit the desired edits in
   chat. Review the sealed design, adaptive thresholds and crystal self-control
   failures. No default rule requires every observed anchor.
6. Send `/guided` to run the **entire configured catalog** for every selected
   template. This has no library pilot prerequisite and no candidate Top-K.
   Independent ChEMBL controls are separate diagnostic molecules, not a sample
   of the library and never included in catalog coverage counts.
7. Use `/status` and `/results`. After completion ask to retain molecules matching
   particular anchors, ALL or ANY, at the original threshold. Selection uses
   one actual pose at a time; molecule-level unions do not create a new pose.

The optional-contact-only design requires at least one adopted anchor in a pose.
Mandatory contacts and each alternative group are additionally enforced in
that same pose. Optional contacts contribute weighted scores. User-reviewed hard
exclusion spheres reject poses; soft spheres penalize scores. These constraints
are separate from receptor heavy-atom clashes. The LLM cannot invent exclusion
coordinates. Supply the coordinate frame, center, radius, evidence and rationale
in chat when such a region is available.

## Funnel and reporting

Every template checks every conformer through: reference-derived heavy-atom
range, typed feature counts, principal extents, anchor necessary conditions,
seed geometry, receptor/hard exclusions, then Gaussian and exact anchor scoring.
Final ranking combines same-pose shape/color Tanimoto, weighted optional scores
and soft exclusions. If no manual pose cutoff is supplied, a baseline is derived
from passing crystal identity-pose scores with a margin and the permissiveness
coefficient. These are exploratory thresholds, not calibrated universal recall.

`report.json` points to each template report, per-stage counts/timings, the union
SQLite database and independent-control report. Each retained pose has its
template, transform, anchor assignments and score breakdown. Exact Murcko
scaffold grouping is used; it is not advertised as 3D similarity clustering.
Increasing template count increases work: this implementation runs templates
sequentially with parallel conformer chunks and repeats catalog I/O. No speedup
or finish time is claimed before measuring the workstation run.

ChEMBL controls require exact target/species, single-protein target, binding assay,
target confidence 9, equality relation and nM Ki/Kd/IC50 at or below 1000 nM.
Assay types and documents remain separate; values are not combined into a single
potency. Reference-cohort chemical graphs are excluded from the independent
panel. Deterministic ETKDG conformers use the library feature-direction function,
but this does not validate the original MOL2 ingestion route. Missing service
data or failed preparation is reported as unvalidated, never as passed recall.

## Boundaries requiring explicit follow-up

Local evidence is recorded in `data/e053-mdm2-implementation.json`: 372 tests
passed, two skipped; 48 admitted instances from 45 PDBs; 42 unique prepared
ligands, 126 contact modes; 43/48 crystal identity controls passed the injected
implementation design. A 24-molecule ChEMBL panel excluded one reference overlap;
the remaining 23 passed the exploratory funnel. This was a deterministic injected
recommendation check, not a live LLM recommendation or full-library validation.

- Transformed biological assemblies, ambiguous interfaces, other receptor states
  and polymer-ligand chemical preparation remain pending; they are not pooled.
- Resolution/occupancy/completeness and geometric admission do not establish
  ligand-density quality or energetic necessity.
- Throughput is measured and the permissiveness coefficient is editable; there
  is no automatic optimizer that sacrifices known-active retention to meet a
  fixed downstream budget.
- This delivery ends at candidate/pose selection. The older WEE1 docking adapter
  is not silently applied to MDM2. General-target SDF export, docking and affinity
  models require their own validated adapters.
- Local fixture tests and real MDM2 structure preparation are not proof of
  25,813,808-conformer coverage. The workstation report must establish that.
