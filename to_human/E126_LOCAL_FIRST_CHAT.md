# E126: conversational interface and local-first scientific audit

Date: 2026-11-10. Implementation is in the isolated E125 delivery checkout,
`D:/agent/work/e125-github-delivery`, based on commit `7dcc8b2`. The local delivery
branch is `delivery/e125-general-20261029`; its GitHub delivery target is
`feature/structure-guided-chat`. The older, dirty development checkout was left
intact. For release installation and acceptance, see [the update guide](E126_UPDATE_AND_TEST.md).

## Completed interface work

1. **Streaming answers and tool progress.** Ordinary natural-language messages
   use `/api/message/stream`. Actual provider content deltas become answer
   previews; tool selection/arguments execute only after complete structured
   validation. Invalid plans or unsupported evidence citations reset the preview.
   Final responses and tool traces remain local and persist across page reloads.
   Closing the browser does not cancel an already submitted scientific task.
2. **Safe Markdown.** Assistant answers support headings, paragraphs, basic lists,
   emphasis, blockquotes, tables, links and fenced code with copy controls.
   Raw HTML is escaped; remote images are not fetched. This is a small supported
   subset, not full nested CommonMark, mathematical typesetting or syntax coloring.
3. **Local file attachments.** Choose or drop structures, sequences, SMILES tables
   and reports onto the composer. Files are preserved byte-for-byte under the
   active Project, with SHA256, original name and session ownership. Maximum:
   eight files/message, 50 MiB/file. Uploads alone do not launch computation.

GPT/DeepSeek profiles use Chat Completions SSE when streaming is requested. For
an endpoint without SSE support, set `"stream": false` in its existing JSON
profile; tool progress remains available, and the final answer arrives together.
There is no automatic provider fallback or resubmission after a transport error.

Attachment bodies are not automatically inserted into the provider prompt. The
model receives local path/name/size/hash references; registered tools determine
what can subsequently be read or processed. Uploading a file is not universal
support for every file format in every scientific workflow. Receptor PDB/mmCIF
and receptor-frame reference SDF/MOL2 paths can enter the generic target search.
An AF3 sequence/SMILES file still needs the existing preparation adapter or the
standalone batch tools; the optimized batch launcher is not a Chat tool yet.
Use natural-language messages with attachments; slash commands use explicit paths.

Keep complete libraries on the Agent host and refer to their configured paths.
An upload from a laptop to a remotely hosted Agent is stored on that Agent host,
not necessarily on the laptop. Input directories outside the Project still require
the existing explicit `--input-root` configuration. Scientific input copying,
task provenance, compute gates and result-import scope are unchanged.

Ordinary chat snapshots no longer import RDKit indirectly just to hash or save
files. Stable polling snapshots also avoid rebuilding the interface every two
seconds. Native scientific dependencies load when the relevant workflow needs them.

## AF3 independent-ligand batches

Each table row remains a separate protein-ligand complex. Explicit homomer or
heteromer assemblies, fixed cofactors, seeds and literal ligand chemistry are
preserved. Matching protein MSA/templates are prepared once and reused; JAX cache
configuration remains in the existing school runner.

The runner now supports explicit restart from a saved report:

```bash
bash tools/af3_smiles_batch/start_batch.sh --resume /absolute/batch/run_previous/report.json
# Equivalent direct invocation:
python tools/af3_smiles_batch/run_batch.py /absolute/batch --run \
  --resume /absolute/batch/run_previous/report.json
```

Resume checks the batch manifest, original/submitted inputs, helper/runner code,
configured runtime paths, final models and available confidence files. Completed
verified jobs are skipped; failed/unrecorded jobs rerun. Shared preparation is
reused only with a matching saved hash and protein identity. Reports are atomic,
attempt directories are unique, and one process lock protects each batch.

This is receipt-based reuse. Interruption before the completed job receipt is
saved can still repeat that job. Runtime path equality does not fingerprint
container, weights or database contents; keep them immutable during recovery.
Legacy reports without the contract require a fresh run. Host checks still do
not count as inference. Predicted complex confidence is not binding affinity.

## 3D screening audit

| Area | Verified implementation | Consequence |
| --- | --- | --- |
| Local coverage | Full catalog iteration, bounded chunk workers, sealed reusable chunk outputs | The existing engine can reuse completed scan chunks and retain results locally. |
| Pose sampling | Finite rigid seeds against the available library conformers | Full catalog coverage is not exhaustive conformational or pose sampling. |
| Early rejection | Size, feature coverage and extent eligibility precede detailed geometric scoring | Some gates are heuristic; reducing survivors does not by itself establish accuracy. |
| Required contacts | Requested contacts are checked within the same candidate pose | Contacts from separate poses cannot be merged into an invented simultaneous match. |
| Retained pool | `candidates.sqlite`, `candidate-poses.jsonl`, cluster/member tables | Keep global, molecule and conformer IDs, source records, masks, transforms and scores. |
| Aggregation | Best representative per molecule, exact anchor mask and spatial-group signature | Aggregated output is not every surviving conformer. Per-chunk pose JSONL preserves the earlier passing records. |
| Scaffold groups | Non-stereochemical Murcko grouping; members retained | Grouping is for organization, not identity equivalence or deletion of stereoisomers. |
| Structure handoff | Budget workflow exports original and transformed MOL2 with provenance/QC | Generic `target_search` currently ends at retained candidates; its automatic MOL2 export adapter is still missing. |

The audit did not change scientific thresholds, candidate quotas or ranking rules.
Existing regression checks cover same-pose logic and source/seal contracts. No
new production recall, enrichment or throughput measurement was made.

## Prioritized local-library work

**P0: connect the generic retained pool to a typed structure handoff.** Accept an
owned completed search task and an explicit selection policy. Reuse verified source
records; export both original-coordinate and receptor-frame structures, chemistry
QC, rejected/review-required records, and a downstream manifest. Preserve the full
pool and pose alternatives. Make molecule versus conformer quotas explicit; do
not silently substitute a fixed Top-K or one conformer per molecule. This closes
the most important current gap between a generic search and the next calculation.

**P1: index source structure locations.** Existing MOL2 export scans each source
sequentially to locate requested records. A versioned byte-offset/length index,
bound to source hashes and molecule/conformer IDs, would allow grouped random
reads across repeated export pages. Build it once per immutable library version;
benchmark the storage device before choosing worker counts.

**P1: make the local candidate catalog the durable interface.** Reuse SQLite and
small manifests rather than repeatedly embedding complete reports in Chat state.
Expose counts, stage timings, checksums, export pages and download/path controls;
load detailed poses on demand. Raw library files and scientific arrays stay in
their existing local store. Separate exact identity, stereochemistry, charge and
protonation policy from scaffold grouping.

**P1: validate narrowing quality and cold/warm performance together.** For each
target, use a recorded known-active/reference set plus an explicitly defined
background set. Compare the deployed funnel with relaxed coarse gates and a more
thorough reference search on an agreed subset. Record molecule/conformer counts,
known-active retention and rejected-case audits, pose recovery where a reference
exists, stage times, preflight hashing/I/O, peak memory, and cache status. This is
a proposed protocol, not an experiment performed in E126. Keep calibration and
held-out evaluation separate. Do not infer recall from an accepted-file count.

**P2: connect shared-MSA AF3 batching to Chat after candidate selection.** Use
explicit verified protein/base inputs and the selected ligand manifest, bounded
GPU scheduling and the new restart receipts. Keep independent complexes separate.
AF3 is a downstream structural examination path, not the present full-library
retrieval engine. The optimized CLI is available now; Chat dispatch remains work.

## Verification and delivery limits

Protocol: `experiments/E126-local-first-chat.md`. Tests use synthetic local inputs,
mock provider SSE and mock AF3 launches. They do not invoke live LLM APIs, AF3 GPU
inference, docking or a production library scan.

197 Python cases passed across the targeted Chat, AF3, generic-query, full-library,
guided, result-import, campaign, planner, structure and budget suites. Four Node
Markdown checks passed; JavaScript syntax and diff whitespace checks passed. These
are targeted checks across runs, not a single full repository test invocation.

Initial broader regressions found two queue tests using root Python without RDKit:
their fixtures replace PYTHONPATH for child processes. Re-running those tests with
the existing scientific virtual environment passed; no production guard was relaxed.
Native chemistry imports were slow on this host, but completed.

The release follow-up provides a working portable temporary preview server at
`scripts/qa_chat_features.py` and an HTTP upload-to-real-tool-to-stream-to-saved-state
integration check. Browser visual acceptance remains pending: browser-extension
access was unavailable and native Edge could not be restored from minimization.
Validate the actual GPT/DeepSeek endpoint, file-picker/drop interaction and viewport
appearance using the update guide after restarting Chat on the workstation.
