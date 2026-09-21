---
name: aidd-3d-search
description: Run and assess this repository's calibrated WEE1 full-library 3D search, Gaussian refinement, E031 annotations and E037 performance comparisons.
---

Use the repository root as the working directory. Read `research-state.yaml` for
current evidence before selecting an experiment or claiming acceptance.

## Choose the execution path

For counts over every library conformer under classified feature conditions, use
the [E044 guide](../../../to_human/E044_FULL_LIBRARY_CONDITIONS.md) and chat
`/full_count CLASSIFICATION_TASK_ID`. Do not substitute exhaustive descriptor
Top-K retrieval: it does not evaluate all poses/features. E044 has no Top-K or
Gaussian Top-N membership budget. It retains one heuristic pose per conformer,
not all orientations/torsions. All/any selection is evaluated within a pose before
molecule deduplication. Pilot/partial results are never full-library counts.

For exhaustive descriptor coverage and classified crystal-derived query features,
read the [E043 guide](../../../to_human/E043_EXHAUSTIVE_CLASSIFIED_SEARCH.md).
New chat searches default to exhaustive retrieval; the low-level CLI needs
`--exhaustive`. Its actual scan belongs to online execution time. Classification
uses native feature hypotheses, requires no PLIP and does not validate candidate
protein contacts. Halogen-specific matching remains unsupported. Preserve original
E031 rankings. Pocket preparation and reference-gated Glide execution are separate
explicit tasks; real chemistry and licensed execution require workstation checks.

- For a single calibrated search, use `scripts/run_e036_fast_3d.sh` or the
  `search_3d` prompt action. Read the [E036 guide](../../../to_human/E036_FAST_3D_RUN.md).
- To compare scheduling and bounded seed generation, use
  `scripts/run_e037_workstation_suite.sh`; read the
  [E037 guide](../../../to_human/E037_ONE_SHOT_VALIDATION.md). Do not also run E036
  separately or repeat the passed E035 million-conformer stress test without a
  new reason.
- For molecule disagreement and pose inspection, use the existing E035 artifacts
  and [review guide](../../../to_human/E035_WORKSTATION_RUN.md).

The prompt adapter accepts only `wee1_qt9`, `wee1_824` or `wee1_both`. The low-level
module `aidd_agent.fast_3d_search` accepts `--query 8bju|1x8b|both`. A new target needs
its own prepared and validated query; never substitute a WEE1 query for it.

## Preserve scientific meaning

Keep the calibrated budget10000, nprobe128, per-objective Top5000 and512pair-seed
cap unless the task explicitly calls for a separately evaluated change. Bounded
generation retains the same unique seed prefix; it does not lower the cap.
E031 is annotation-only. Molecule caps and budget retention are not chemical
rejection, activity enrichment or pose-quality acceptance.

Validate query/library identities and source hashes through the adapters. Use new
output directories when code, parameters or inputs change. Completed receipts
preserve original timing; partial reuse and profiling are ineligible for fresh
latency comparisons. Integrity scans/index loading and exact-reference scans are
not online query latency; summed worker compute is not wall time.

Return the actual report paths, gate states, measured timing scope and any failed
equivalence checks. Distinguish user-reported workstation evidence from local
fixtures. Inspect exported discordant poses before interpreting ranking differences.

## Human selection and docking handoff

Use the [E042 guide](../../../to_human/E042_SCREENING_TO_DOCKING_HANDOFF.md) for
chat evidence review, explicit same-pose selection previews and separately confirmed
pose/ID exports. These are session-bound local coordinator actions, not actions
the model may append to a search plan. Reuse completed search artifacts and their
receipts. Never invent anchor IDs or default a selection threshold. HBA/HBD feature
matches are hypotheses, not validated candidate hydrogen bonds. A user selection
branch changes membership only after explicit criteria; original E031 ranks stay
unchanged. Export is a preparation handoff, not a docking result.
