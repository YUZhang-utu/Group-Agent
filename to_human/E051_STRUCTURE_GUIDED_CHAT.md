# Structure-guided screening through chat

This is a new, isolated workflow. Keep the running E050 checkout and outputs
unchanged. Use a second checkout and a separate chat storage root/port. Do not
start competing full-library jobs while the existing task occupies the machine.

## Install the branch without changing the running checkout

Run these once on the workstation:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git fetch origin feature/structure-guided-chat
git worktree add --detach ../Group-Agent-e051 origin/feature/structure-guided-chat
cd ../Group-Agent-e051
```

Use the existing scientific Python environment, including RDKit, Biopython,
NumPy and SciPy. Check dependencies without changing that environment:

```bash
python -c 'import rdkit, Bio, numpy, scipy; print("Scientific dependencies available")'
```

Reuse the existing runtime profile (its `search` section specifies the library
batch and worker counts) and LLM profile directory. Set `AIDD_RUNTIME_PROFILE`
and `AIDD_LLM_CONFIG_DIR` to those actual paths if they are not already exported.
Keep API keys in their existing environment variables.

```bash
bash scripts/run_chat_agent.sh \
  --storage-root /mnt/local/hand/yuzhang/aidd/e051-chat-workspace \
  --runtime "${AIDD_RUNTIME_PROFILE:?Set the existing runtime profile path}" \
  --llm-config-dir "${AIDD_LLM_CONFIG_DIR:?Set the existing LLM profile directory}" \
  --port 8766 --allow-compute
```

Open the chat at `http://localhost:8766` on the workstation, or use your existing
SSH forwarding method. This command starts the chat service; it does not launch
a full-library scan by itself.

## Conversation sequence

1. Ask: "Survey ligand-bound PDB structures for human WEE1. Verify the target,
   analyze crystal interactions and recurring contacts, and show the evidence.
   Do not start a library scan yet." Other targets are supported through verified
   UniProt identity; ambiguous gene/organism identity requires clarification.
2. After the survey finishes, send `/recommend`. The configured GPT or DeepSeek
   provider receives bounded structural evidence and proposes a small mandatory
   core, alternative groups and optional anchors, with citations and uncertainty.
3. Review/edit the proposal in natural language using the displayed anchor IDs,
   for example "Make anchor [ID] mandatory, make [ID] optional, and keep the
   remaining recommendations." Or send `/adopt` to use the proposal as written.
   An edit creates an adopted design; it does not start screening. The original
   crystal ligand must pass the anchor and pocket controls before adoption.
4. Once machine resources are available, send `/guided`. This scans the entire
   configured library; there is no pilot prerequisite or Top-K/Top-N cap.
5. Use `/status` or `/results`. Each completed task posts its next action.
   Commands accept an optional task ID when several designs exist.
6. After screening, ask "Keep molecules with both [anchor ID] and [anchor ID],
   at the original score threshold 0.5." The selection runs against real saved
   poses and produces molecule IDs and pose JSONL. It preserves scaffold IDs.
   Anchors found on different poses never count as a simultaneous match.

Without a usable crystal, the survey searches Europe PMC and distinguishes API
failure from no hits. `/recommend` summarizes the available literature. Provide
a target-bound PDB ID in chat for a new survey. Arbitrary uploaded receptor/ligand
files are not yet a chat ingestion adapter; literature cannot create coordinates.

## What the new funnel does

Full conformer catalog -> heavy-atom size -> typed feature coverage -> shape
extent -> anchor necessary conditions -> mandatory feature geometry -> original
seed generation -> same-pose geometry feasibility -> pocket exclusion -> Gaussian
ranking on surviving seeds -> exact anchor assignment -> molecule aggregation
and exact Murcko scaffold groups -> human anchor-combination selection.

Mandatory anchors must all match; each alternative group requires at least one
member; optional anchors are reported without becoming rejection criteria.
All requirements apply to one pose. Removing an early mandatory rule requires a
new design and scan: discarded molecules cannot be recovered by post-selection.

The pocket default rejects any ligand heavy atom closer than 1.2 Angstrom to a
protein heavy atom. This catches severe overlaps only. It is not van der Waals
validation, water/metal treatment, pocket flexibility or docking. Thresholds are
versioned exploratory defaults, not calibrated target-specific recall. Failure
of the crystal reference blocks adoption and reports preparation needs.

Repeated contacts are mapped through target-associated protein entities to the
verified UniProt sequence. Counts use distinct PDB entries and distinct CCD
ligands; identical copies do not increase recurrence. The default survey analyzes
at most 12 structures; all discovered and unanalyzed entries are reported.
Ambiguous mapping and ligand alternate locations are reported, not pooled.
Crystal-feature matches remain interaction hypotheses, not binding proof.

## Reports and verification

The chat task's report path leads to the survey, recommendation or adopted
design. A guided task includes `guided/full-library/suite-report.json`, progress,
sealed chunks, `candidate-poses.jsonl`, `cluster-members.csv`, `clusters.csv` and
`candidates.sqlite`. Post-selection creates `selected-poses.jsonl` and
`selected-molecules.txt`. The new selection is not yet an SDF/docking adapter.

Inspect `guided_geometry_passed/rejected`, `original_seeds`, `possible_seeds`,
`pocket_tested_seeds`, `pocket_rejected_seeds`, `pocket_passed_conformers`,
`gaussian_evaluated_seeds`, `gaussian_evaluated_conformers`, final matches and
stage timings. Guided substage counts are conformer/seed counts; the existing
eight-level molecule-stage table groups the geometry/pocket gates together.

This version changes Gaussian ranking scope to surviving seeds. Gaussian remains
a ranking operation, so equal Gaussian/matching conformer counts can still be
legitimate. Speed and retention improvements must be measured on the workstation;
no specific rejection percentage or speedup is promised. Local fixture tests do
not establish real full-library recall or live API/LLM reliability.
