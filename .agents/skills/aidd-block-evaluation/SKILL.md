---
name: aidd-block-evaluation
description: Sample frozen macrocycle partitions, prepare shared PLANTS jobs and analyze block-aware docking results through the existing prompt agent.
---

Read [the E097 protocol and workstation guide](../../../to_human/E097_BLOCK_PLANTS_EVALUATION.md).
Use the configured E094/E095/E096 artifacts, preserving source conformer IDs and all
scheme/block memberships. Sample conformers uniformly without replacement; deduplicate
docking inputs by conformer, not molecule. Preserve raw coordinates when exporting.

For PDB-only requests use protein_from_pdb and receptor_assess to identify the target
and compute the existing evidence-based receptor recommendation. Wait for explicit
human acceptance before adoption. The adopted structures can be cleaned and converted
with plants_receptors using the configured SPORES executable and its installed mode.
Compute the native ligand centroid/radius and preserve the receptor coordinate frame.

Use block_sample, block_plants_prepare, block_plants_run and block_analyze with typed
preceding-step references. Existing-task follow-ups use the owned block_evaluation
chat coordinator. Trusted runtime profiles provide paths, prepared receptors and
binding sites. Do not infer receptor coordinates from a target name. Preparation
does not authorize execution; respect the user's requested stages and existing compute
authorization. Missing receptor inputs do not prevent sampling/export.

Report actual receipts and failures. Local synthetic fixtures do not validate live
PLANTS, MDM2 enrichment, binding affinity or retrieval recall. Rigid ligand mode is
the default conformer comparison; flexible docking changes its interpretation. Preserve
pose identity for N-E rescoring. Never reject unsampled blocks from this pilot alone.
