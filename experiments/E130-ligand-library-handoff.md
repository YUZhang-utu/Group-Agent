# E130: ligand-led whole-library search and MOL2 handoff

Date: 2026-11-11. User clarified two entries: target-driven experimental complex
discovery/ligand analysis, or an attached small-molecule structure, followed by
3D search against the complete local library and candidate MOL2 delivery.

Protocol: preserve the existing target discovery/diversity/consensus and contact-
first budget workflow. Add a trusted ligand-library tool for session uploads or
explicit prepared-query IDs from owned completed structure tasks. Reuse Gaussian
chunk scoring over all catalog conformers, deduplicate by molecule with stable
best-score/conformer ties, and reuse original-record MOL2 validation/export.
Do not invent receptor contacts for a standalone ligand. Store complete ranking
locally; delivery count is an output cap, not a retrieval cap. Seal inputs/config
and retain resumable chunk receipts. Test upload ownership, source-query selection,
full scoring coverage, identity/transform preservation, local export and failure
propagation on a small synthetic library. No production speed/recall claim.

## Results

93 cases passed in the integrated search/Chat/planner/budget/chain suite. The
small-library fixture builds actual conformer and chemical artifacts and uses
real Gaussian scoring and original MOL2 export, including a real Chat worker
subprocess for the uploaded path. Prepared cocrystal discovery is represented by
a saved-query fixture, not live PDB evidence. Six focused cases cover both entry
paths, best-conformer deduplication, ownership, compute gating, 2D rejection,
queue receipts, source changes, catalog identity, chunk reuse and output paths.

An additional lineage check initially assumed a chemical-catalog top-level
conformer count; the actual format records per-shard ranges. Corrected to compare
range coverage and the source artifact-catalog checksum. Final tests passed with
the real catalog builder. Four Markdown checks and JS syntax passed; no browser
visual check or production-library timing/recall run was performed.

Final affected subset: nine tests passed, including three additional workflow
navigation cases (96 distinct passing cases across runs). English guard: 743 files.
