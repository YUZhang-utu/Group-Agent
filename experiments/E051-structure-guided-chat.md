# E051: structure-guided chat workflow

Status: implementation in isolated feature/structure-guided-chat worktree.
The running E050 checkout, outputs and full-library task must not be modified.

User request: chat-driven protein/ligand PDB discovery, interaction evidence and
recurrence, evidence-grounded LLM advice with user edits or explicit adoption,
literature fallback when crystals are missing, and stronger pre-Gaussian
geometry/pocket filtering feeding the full-library preselection workflow.

Workflow: verified target -> structure survey -> advisory proposal -> human edit
or adoption -> immutable search design -> guided full-library scan -> pose and
scaffold tables. Silence is not an approval event. A user may explicitly delegate
design choice to the LLM; the default proposal can then be adopted in chat.

Survey records all discovered entries and bounded analyzed structures separately.
Repeated contacts count distinct structures, not crystallographic copies. Residue
recurrence needs target sequence mapping; unverified mappings are not pooled.
Direct hydrogen bonds without explicit hydrogen angles and other geometric
contacts are hypotheses, not experimentally proven energetic necessities.
Literature records retain citations/abstract provenance but cannot manufacture
coordinates. Missing/failed APIs are not proof of no structures or no papers.

Geometry uses necessary pair constraints for mandatory anchors and explicit
alternative groups. Pocket filtering applies only to known receptor-frame seed
poses before Gaussian. The reference crystal pose must pass the same pocket
gate; a failure blocks the design instead of silently weakening thresholds.
Novel thresholds are versioned exploratory defaults, not calibrated sensitivity
guarantees. Exact assignment rechecks the adopted rule; optional features remain
reported and do not silently become requirements. Protein clashes are not docking.

Verification: identity/mapping ambiguity, repeated-copy deduplication, no-crystal
fallback and network failure, fabricated LLM IDs/citations, manual overrides,
same-pose rule semantics, pocket transform conventions, reference clash control,
chat session ownership, plan receipts, source seals and full-library adapter tests.
Live RCSB/Europe PMC/LLM and workstation performance remain separate evidence.
