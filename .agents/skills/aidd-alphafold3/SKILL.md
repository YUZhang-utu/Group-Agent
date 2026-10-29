---
name: aidd-alphafold3
description: Prepare and execute local AlphaFold 3 tasks through this repository's verified protein/CCD adapters, then import prediction artifacts, confidence and provenance.
---

Read the [AF3 workstation guide](../../../to_human/E038_PROMPT_PROTEIN_AF3.md).
The implementation is `src/aidd_agent/prediction.py` and the `af3_prepare`/
`af3_run` adapters in `prompt_workflow.py`; reuse them rather than constructing
unreviewed shell commands.

## Prepare

Use literal `protein_sequences` from the user's prompt, a preceding verified
`protein_fetch`/`protein_resolve` result through `protein_step`, or multiple verified
steps through `protein_steps`. Accept an
explicit full sequence or user-specified construct interval; the LLM must not
generate the sequence. The current prompt adapter supports a monomer or explicitly
requested homomer (`protein_copies`, integer 1..8, default 1), or heteromeric sequence
lists with repeated entries for explicit stoichiometry, and
optional RCSB-validated CCD ligands or user-supplied SMILES ligands. SMILES must be
copied literally and pass strict RDKit parsing; preserve stereochemistry and charge.
Do not infer missing stereochemistry or enumerate protonation states silently.
Basic AF3 dialect-v1 JSON preserves configured
seeds and construct provenance. Installed AF3 handles MSA/templates.
Explicit seed lists support 1..100 distinct integers. `seed_count` optionally
with `seed_start` expands a deterministic range (default start 1); never combine
the list and count modes. Literal sequences are not database-verified identities.

All ligand entries in one input are co-ligands in one complex. Independent SMILES
predictions need separate prepare/run pairs; the shared-MSA batch runner remains
CLI-only. Read-only external batch/manifest import uses `workflow_results` and
does not verify model payloads or create resumable Chat execution tasks.

DNA/RNA, custom CCD/covalent bonds and automatic domain
selection are not exposed by this prompt adapter. Do not silently approximate them
with another task. State the missing capability or prepare a separate implementation
when that is the user's request.

## Execute and import

The trusted local runtime selects native Python or Apptainer execution. Native
profiles use host Python/runner paths; container profiles use host image/model/data
paths and container Python/runner paths. See the
[container guide](../../../to_human/AF3_CONTAINER_AND_LLM_PROVIDERS.md). Reuse the existing
installation and profile checks; do not invent paths, reinstall AF3 or download
weights merely to prepare a plan. Credentials remain outside model output.

Use the existing compute capability when the user has authorized execution;
`--allow-compute` enables the adapter. Without it, preserve prepared inputs and
report `blocked`, not successful inference. Execution uses an argument list with
`shell=False`. Preserve separate failed-attempt directories and execution logs.

Import the final model, structure hash, input provenance and confidence using
`inspect_alphafold3_output`. Multiple/missing final models are failures, not a
reason to pick an arbitrary file. Report model version as recorded, or explicitly
unspecified. AF3 confidence is not binding affinity, docking quality or biological
acceptance. Resume only matching plans/configuration/code and valid output receipts.

Chat confidence explanations and desktop PyMOL controls are documented in
[the structure review guide](../../../to_human/E066_CHAT_STRUCTURE_REVIEW.md).
Local adapter tests do not establish successful AF3 inference or desktop GUI operation.
