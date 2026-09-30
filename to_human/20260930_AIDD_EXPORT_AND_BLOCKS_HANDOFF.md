# AIDD checkpoint: exported MDM2 conformers and running E094

## Resume instruction

Read this handoff and research-state.yaml. Preserve two independent workstreams:
MDM2 conformer selection and E094 work-block construction. The user reports the
3.8-million-conformer export finished and the block job is still running. Do not
restart either operation. No workstation process or final export receipt has been
directly inspected from this host. This supersedes the September 29 running-E091
checkpoint; E091 already finished, and E094 is the currently reported running job.

## MDM2 export and corrected selection unit

User-pasted reports establish complete searches for human MDM2 Q00987:

- PROMPT-59addb9ffd6444c3: 11 templates, 1,000,000 ranked molecules.
- PROMPT-198eeebebde54547: eight templates, 1,000,000 ranked molecules.
- Each exported ranks 1-100,000 as representative conformers, with no shortfall.
- Both retrieved through the complete index (8,318,351 molecules; 25,813,808
  conformers) using USRCAT FAISS ANN. They did not score every library pose.
  Current-target ANN recall is unmeasured.

For the 11-template run, the user confirmed selected-molecules.npy has 1,000,000
entries and selected-conformers.npy has 3,854,146 entries. The latter expands
selected molecules to all stored conformers. The user clarified the desired final
budget is approximately ONE MILLION CONFORMERS, not one million molecules.
They explicitly chose to finish exporting all 3,854,146 first, then select later.
Latest message reports the export completed. Exact exported counts/hashes still
need verification from its report.json; do not call them independently validated.

Search directory:

```text
/mnt/local/hand/yuzhang/aidd/e054-chat-20260923-104806/users/workstation/projects/prj-42c81b1322ab-prompt-aidd/runs/PROMPT-59addb9ffd6444c3/execution/guided/guided/search
```

Original batch: /mnt/local/hand/yuzhang/aidd/library-precompute-20260911.

Expected export directory from the supplied command:

```text
/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-1m-all-conformers
```

Standalone exporter published at commit 37e958b on feature/structure-guided-chat:
scripts/export_selected_conformers.py. Downloaded separately so E094 checkout
need not change. Three local tests passed; source hashes, identity and counts are
checked. Guide: MDM2_ALL_SELECTED_CONFORMERS.md. Outputs: molecules.csv,
manifest.sqlite, part-*.mol2, matching part CSVs, part receipts, report.json.
Coordinates are original library conformers, not docking poses.

Next steps:

1. Inspect export report.json, confirm 1,000,000 molecules and 3,854,146 conformers.
2. Inspect the existing search/ranking.sqlite schema and saved scoring artifacts
   to establish which individual conformers retain scores. Do not assume every
   conformer score survives merely because one million molecules were ranked.
3. Design selection of approximately one million unique conformer IDs from the
   exported set; multiple conformers per molecule are allowed. Preserve rank,
   template, score scope and source mapping. If scores are only representatives,
   state that limitation and propose an explicit policy before claiming a global
   top-million conformer ranking.
4. Do not take the first million selected-conformers.npy entries: these are global
   IDs, not score order. Do not rerun ANN or regenerate the original MOL2 export.

## E094 block job

User reports it is running, with an earlier rough nine-hour estimate. There is no
measured finish time or direct process access. Keep code, parameters and outputs
unchanged while it runs. Published build code: 79b9d4e.

Build input:
/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone

Expected output:
/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/final-blocks-e094

Scope: 24,663,736 admitted conformers. Regular execution blocks have a 5,000
minimum and no capacity ceiling; chemical identities remain traceable within
pools. Side-chain property/steric refinement requires full profile coverage,
minimum child size 5,000 and at most eight children per parent, without forced
splitting. E091 special pool remains 132,885 (0.538787%), pending search evaluation.
The separate 1,150,072 source chemistry-review records remain deferred.

After completion, inspect backbone-validation.json (structural_gate passed),
properties-validation.json (structural_gate and property_gate passed), and
properties/report.json (blocks_below_minimum zero and conserved counts). Guide:
E094_FINAL_WORK_BLOCKS.md. Passing establishes offline work-block integrity,
not search recall or docking usefulness. Later evaluation uses 3D search,
PLANTS/MDM2 and the user's N-E models; affinity model is still undecided.
