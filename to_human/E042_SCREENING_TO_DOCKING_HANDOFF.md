# Prompt-driven screening evidence and human selection

Recovery from `Retrieval/refinement molecule identity mismatch`: update to the
bytes/text ID compatibility fix and restart the chat server. Retrieval stores S16
byte IDs while Gaussian refinement stores U16 text IDs. The review now decodes
IDs before comparing them; genuine mismatches still fail. Keep the completed
search. Submit a NEW `/evidence SEARCH_TASK_ID` using that search's chat task ID,
not the failed review task ID. Do not resume the failed review after a code update,
because its execution protocol is bound to the previous code. Do not edit NPZs or
receipts to bypass validation.

This workflow has four chat turns: search, evidence review, selection preview and
explicit export. Reuse a completed search when possible. No docking calculation
is submitted by these steps. The implementation supports the calibrated WEE1
queries; another target still requires query preparation and calibration.

## What the existing search means

USRCAT retrieves a fixed descriptor-ranked budget from the full indexed library.
Gaussian coarse selection allocates refinement to a union of objective Top-N
candidates. The anchored Gaussian objective weights selected crystal-ligand
features. E031 subsequently annotates stored poses using typed, directional,
one-to-one feature assignments. E031 does not change the original admission or
ranking. Budget reduction is not a measured chemical rejection rate.

The original anchor extractor uses ligand donor/acceptor features, complementary
protein atom roles and a maximum 3.5 A feature-to-partner distance. It records the
nearest compatible partner and burial estimates. The baseline extraction does not
evaluate explicit-hydrogen angles. Read each actual query manifest: later query
preparation may contain additional direction/angle evidence. An anchor is a
ligand feature with recorded crystal-contact evidence, not a protein residue or
a verified new hydrogen bond in every candidate. Candidate anchor scores measure
feature reproduction, not binding probabilities. Receptor clashes and flexible
pose quality are separate questions.

The review groups anchors by ligand HBA/HBD feature class and includes the crystal
query ID, zero-based ligand atom indices, protein chain/residue/atom, distance,
angle status, burial and original evidence fields. Hydrophobic, pi-stacking,
salt-bridge, water-bridge and metal-coordination classes are explicitly unassessed.
Legacy query mappings are reconstructed from the small crystal query and accepted
only when every stored Gaussian query array is identical. This does not rescan
the library or repeat Gaussian candidate refinement.

## Start on the workstation

```bash
git pull --ff-only origin main
export AIDD_LLM_CONFIG_DIR=/mnt/local/hand/yuzhang/aidd/config/llm
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime-apptainer.local.json
bash scripts/run_chat_agent.sh --allow-compute
```

Keep the selected provider's API key in its existing environment variable. Open
the printed loopback URL. Stop an older chat server before starting this version.
The UI, reports and prompts below are English. The model may interpret other input
languages, but replies remain English.

## 1. Search or attach the existing completed search

For a fresh search, enter:

> Run the calibrated full-library WEE1 search using both QT9 and 824 queries. Keep E031 annotation-only. Do not run docking or choose selection thresholds.

Alternatively use the existing-plan attachment field with the `plan.json` next to
the completed run's execution directory. For the previously reported run:

```text
/mnt/local/hand/yuzhang/aidd/prompt-workspace/users/workstation/projects/prj-252fa94197d6-prompt-aidd/runs/PROMPT-4689d49ad1834f89/plan.json
```

Attachment does not rerun the old code. Use the task ID displayed by the chat UI,
not the PROMPT directory name, in the following requests.

## 2. Inspect evidence and counts

> Review the completed search task TASK_ID. Show library, retrieved and refined conformer/molecule counts. Classify its crystal-derived anchors and show each anchor's ligand atoms, protein partner and geometry evidence. Do not select or export candidates yet.

Deterministic equivalent: `/evidence TASK_ID`. After completion ask:

> Where are the evidence table and results for this task?

Or `/results REVIEW_TASK_ID`. The report lives under
`execution/screening/screening/report.json`; alongside it are `anchors.csv` and
`8bju-candidate-anchors.csv` / `1x8b-candidate-anchors.csv` for queries present.
The latter contain per-conformer scores for the stored anchored-objective pose.
Counts at 0.25, 0.5 and 0.75 are diagnostic examples, not calibrated acceptance
thresholds. Per-query counts overlap: the report also provides molecule unions.
Molecules use existing source-grouped identity, not newly deduplicated chemical
structures across vendors.

## 3. Request an explicit selection preview

Copy actual full anchor IDs from the evidence response; do not use placeholder IDs.
For example, after choosing two anchors and a threshold yourself:

> From review task REVIEW_TASK_ID, preview molecules matching all of these anchors: EXACT_ANCHOR_ID_1 and EXACT_ANCHOR_ID_2. Use minimum anchor score 0.5. Do not apply a molecule cap. Do not export or dock yet.

Here 0.5 is an illustrative user-selected setting, not a scientific recommendation.
Use `any` instead of `all` for an OR condition. If a cap is wanted, specify it
explicitly, for example `at most 100 molecules`. Every preview uses anchors from
one crystal query. Make a separate preview for the other query; cross-query
intersection/union selection is not implemented. Never combine coordinates or
anchor matches from different crystal frames.

ALL conditions must hold in the same stored conformer/pose. After qualifying,
deduplicate molecules and keep the highest original Gaussian-scoring qualifying
pose, breaking ties by global ID. Apply an optional molecule cap last. The report
distinguishes matching conformers, matching molecules, selected molecules and
selected representatives; `selected-molecules.csv` records IDs. Original search
files and rankings remain unchanged. A zero-match result is valid and is not
silently relaxed. This is an explicit human selection branch, not an automatic
E031 rejection gate.

## 4. Confirm the preview and export

> Confirm selection task SELECTION_TASK_ID and export its selected molecule IDs and representative poses for docking preparation. Do not submit docking jobs.

Deterministic equivalent: `/export SELECTION_TASK_ID`. This must be a separate
turn after the preview. The export contains:

- `selected-poses.sdf`: one qualifying rigid pose per selected molecule.
- `selected-ids.json`: source-grouped molecule, conformer, global ID and query.
- `report.json`: counts, exact policy, source/output hashes and preparation status.

Ask `/results EXPORT_TASK_ID` for paths. The SDF contains heavy atoms, bonds and
formal charges in the original crystal frame. Full stereochemical annotation is
not reconstructed. It is not automatically a prepared Glide or PLANTS input.
Ligand chemistry/protonation, hydrogens, receptor preparation, site/grid and a
reference-ligand docking validation are still required. Installed executables
alone do not establish these steps. No binding activity claim follows from export.

## Acceptance checklist

1. Library counts match accepted evidence, and retrieval/refinement counts match
   the actual saved arrays. Independent queries are not summed as unique molecules.
2. Manually inspect at least one HBA/HBD anchor and its named crystal partner.
3. Try ALL versus ANY on the same anchors. ALL must not exceed ANY before caps.
4. Try a strict threshold; preserve zero results rather than inventing candidates.
5. Verify preview count equals the number of exported IDs and SDF records.
6. Confirm status/result questions and evidence/selection tasks did not rerun search.
7. Changing a source artifact must fail the dependent task; another session's task
   ID cannot be selected or exported.

Local synthetic fixtures validate these engineering contracts, including a real
background subprocess chain. Live DeepSeek/GPT routing, real library exports and
manual pose inspection remain workstation acceptance tasks. Compact evidence
summaries (anchor IDs, residue geometry and counts) are sent to the selected LLM
to interpret follow-up selections. Raw library structures, arrays and logs are
not sent by this feature.
