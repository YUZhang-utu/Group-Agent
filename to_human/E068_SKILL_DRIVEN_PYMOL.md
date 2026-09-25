# Skill-driven desktop PyMOL in Chat

The existing GPT/DeepSeek profile now plans composable PyMOL API calls using the
Google DeepMind science-skills PyMOL skill and live task-object metadata.
This extends the desktop bridge; no new model service or PyMOL installation is needed.

Upstream: https://github.com/google-deepmind/science-skills/tree/68832757cbbf941c620b71df5756cf6e5cc287b0/skills/pymol

The unmodified skill, two reference documents, Apache-2.0 license and content hashes
are bundled under `src/aidd_agent/skill_data/pymol/`. Hashes are checked when planning.
Please check the license applicable to your PyMOL installation at https://www.pymol.org/.
The upstream headless setup is adapted to reuse the persistent desktop process:
the model cannot launch, install, load arbitrary files or quit PyMOL.

## Update and start

1. In the workstation repository, run `git pull --ff-only origin feature/structure-guided-chat`.
2. Stop the existing Chat server and close its bridge-controlled PyMOL window.
   Save any manual work first. Both processes must restart to use the new code.
3. If the project is installed editable, no reinstall is required. Otherwise install
   the updated package using the same environment as the previous Chat service.
4. Reuse the existing launch command and environment variables:

```bash
bash scripts/run_chat_agent.sh \
  --storage-root "$AIDD_NEW_WORKSPACE" \
  --runtime "$AIDD_RUNTIME_PROFILE" \
  --llm-config-dir "$AIDD_LLM_CONFIG_DIR" \
  --port 8766 --allow-compute
```

5. Open the existing conversation and choose GPT or DeepSeek with working credentials.
   Open a completed structure task; no scientific search needs to be repeated:

```text
/pymol YOUR_COMPLETED_TASK_ID
```

Wait for the desktop structures to load. `/pymol_status` reports the actual receipt.

## Acceptance prompts

The explicit command below bypasses general intent routing. Natural language alone
also routes to the new planner once a task is open.

```text
/pymol_agent Show each protein as a cartoon. Use different carbon colors for the ligands, blue for nitrogen and red for oxygen. Show complete protein residues within 5 angstroms of each ligand as sticks. Label their residue names and numbers. Draw possible polar contacts separately within each complex. Use a white background.
```

```text
/pymol_agent Inspect the protein chains in each complex. Remove only redundant protein chains that are not within 5 angstroms of that complex's ligand. Preserve the ligand. If multiple chains contact the ligand, ask me which to retain.
```

```text
/pymol_agent Color the protein from blue to red along the sequence, keep ligand atoms colored by element, turn the view 30 degrees around y, and zoom into the ligand pocket.
```

```text
/view {"operation":"undo"}
```

Check the desktop appearance yourself. A successful program creates a PNG, PSE,
pre-edit PSE checkpoint, readable API program and execution receipt under
`STORAGE/chat/viewers/SESSION_ID/`. Download buttons appear on the viewer card.
The `*-agent.json` audit records the requested plan, actual receipt, model metadata
and pinned skill provenance. Provider credentials are not written there.

## Execution contract and limits

- The model receives object/ligand identifiers, chain inventories, bounded residue
  examples and API errors. Coordinates and screenshots are not sent to the model.
- Programs are parsed and interpreted as literal calls to an allowed subset of
  `cmd`: representations, colors, selections, labels, distances, alignment, camera
  changes and display settings. No Python eval/exec, imports or shell commands.
- Atom selections are intersected with loaded task objects. `ai_` is reserved for
  generated selections, distance objects and custom colors; do not use that namespace
  for manual objects. Camera/background changes affect the entire viewer.
- Empty selections fail visibly, except count queries which may return zero.
  Execution errors trigger session rollback. The model may repair once after a
  successful rollback. Pending operations are never blindly submitted again.
- Undo restores the session before the last successful agent program, including
  camera state; subsequent manual edits are also reverted. Undo is process-local
  and cleared when opening task structures. The saved checkpoint remains available.
- Polar distance candidates are not a complete interaction analysis. Salt bridges,
  pi interactions and water-mediated contacts need separately validated methods.
- Unsupported requests produce clarification or explicit failure. This is not
  unrestricted PyMOL Python. A new primitive can be added to the API contract
  without adding a new natural-language intent for every display recipe.
- Local fixture tests do not establish real-provider or desktop GUI success.
  Verify these prompts on the workstation before a live demonstration.
