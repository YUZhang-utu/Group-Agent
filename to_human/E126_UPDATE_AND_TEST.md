# E126 update and interface acceptance

## Pull into the existing workstation checkout

Use the existing AIDD Python environment and repository. Stop the Chat server
cleanly before replacing its code, when no scientific job is actively executing.
Do not replace the existing storage root, provider profiles or runtime files.

```bash
git status --short
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
python -m pip install -e . --no-deps
```

If local changes prevent the switch/pull, preserve and reconcile them rather than
resetting the checkout. The update includes the browser assets as package data.
Restart using your existing command and configuration, for example the existing
`scripts/run_chat_agent.sh` launcher. Keep `--storage-root`, `--runtime`,
`--llm-config-dir` and any `--input-root`/`--result-import-root` arguments unchanged.
Keep the existing compute setting; no new GPU job is needed to test the interface.

Open the local access URL printed by the restarted server, including its `#token`
fragment. If viewing through the existing SSH tunnel, retain that tunnel and the
corresponding local port. Refresh the page after restarting; refreshing alone does
not reload the Python server. Select the existing GPT or DeepSeek provider.

## Credential-free interface fixture

For an isolated test that cannot submit scientific jobs:

```bash
python scripts/qa_chat_features.py
```

Open its printed URL. It creates a fresh temporary Project, uses a synthetic model
and calls the real read-only local task-list tool. It does not read your API keys,
real conversation storage or library. It prints the temporary storage location.

1. Send `Inspect my local tasks.` Check tool progress, incremental answer text,
   a formatted table, a bullet list and a code block with a copy button.
2. Attach a small synthetic CSV using **Attach files**, then send a prompt.
   The response must say the attachment references reached the planner. Confirm
   the file is saved under the printed temporary storage with `attachment.json`.
3. Repeat with drag-and-drop onto the composer. Remove an attachment before
   sending to test the attachment chip control.
4. Reload. The final answer and conversation must remain; no duplicate task is
   created. Inspect the Results drawer for the completed research-agent trace.
5. Stop this fixture with Ctrl+C and return to the real Chat server. The fixture
   uses separate storage and must not be mistaken for the configured live Agent.

## Live Agent acceptance

The interface fixes do not replace a working provider/runtime setup. After the
isolated check, use the production server with a read-only prompt first:

```text
Inspect the configured local library status and list this conversation's tasks.
Do not start a search or prediction. Show a concise table of what is available.
```

Confirm the tool trace includes real `library_status`/`tasks` evidence and that no
new scientific job was queued. This confirms the selected provider can select
registered tools and interpret their results. It does not validate GPU inference.
If an OpenAI-compatible endpoint rejects streaming, set `"stream": false` in that
provider's existing local JSON profile and restart. Do not put API keys in Git.

For a molecular-input check, attach your prepared receptor PDB/mmCIF and posed
reference SDF/MOL2. Supply their correct receptor chain identities and explicitly
request **query preparation only**. Inspect the completed query report before
authorizing a full-library calculation. Receptor and ligand must share a coordinate
frame. This does not automatically validate every uploaded format for every tool.

The generic retained-pool-to-MOL2 adapter and optimized shared-MSA AF3 Chat launcher
remain pending; the existing budget export and standalone AF3 batch CLI are
available. Pulling this release does not turn those pending adapters into UI
capabilities. AF3 recovery is documented in `tools/af3_smiles_batch/README.md`.

## Release checks performed here

Previous E126 validation passed 197 targeted Python cases and four Markdown tests.
The release adds an HTTP integration regression covering served JS/CSS assets,
session creation, binary upload, a natural-language request, the real local tasks
tool, streamed answer events, saved conversation reload and an accessible tool
trace. It uses a synthetic provider; no live-provider/GPU acceptance is implied.
The final affected release suite passed all 58 cases in 14.54 seconds; the four
Markdown checks and JavaScript syntax checks also passed.

The temporary interface server starts successfully. Browser automation remains
unavailable: browser access host is absent and the native Edge window could not
be restored from its minimized state. Actual visual/file-picker/drag interaction
acceptance is therefore still to be completed using the steps above. Do not record
those browser interactions as passed based only on HTTP tests.
