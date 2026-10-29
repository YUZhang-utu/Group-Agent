---
name: aidd-target-search
description: Prepare explicit target-independent receptor/reference queries and dispatch the existing full-library guided engine through Chat.
---

Read [E125](../../../to_human/E125_GENERAL_PROMPT_WORKFLOWS.md). Use the existing
`target_search` coordinator rather than inventing shell commands or fixed-target
templates. `prepare` creates a query; `start` prepares and scans; `run` consumes an
owned completed query task. Progress/result questions never submit another scan.

Require explicit receptor chains and a receptor-frame reference ligand, or exact
PDB/CCD/author chain/residue identity. Import exact user-supplied file paths only
from Project or configured input roots. Preserve cofactors explicitly identified
by the user. A supplied receptor is a structural hypothesis, not automatic cohort
or biological validation. Missing site information requires clarification.

Use reported contact classes and residue identities. Per-residue alternatives
and every required residue must be satisfied in one pose. Native controls block
inconsistent receptor/query/policy combinations; do not silently relax them.
The scan covers the full catalog with finite rigid seeds, not all conformations.
Broad coarse defaults, Gaussian ranking and gross overlap checks are exploratory.
Outputs are candidate evidence; docking, affinity and arbitrary-target recall
remain separately validated. Existing sealed run/code/configuration checks apply.
