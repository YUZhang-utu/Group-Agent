# PDB pocket-state prerequisite for new target projects

Status: user-approved workflow design and planner guidance. Pocket clustering,
persisted user adoption and an executor-enforced prerequisite remain to implement.
Existing ligand-diversity and consensus adapters do not implement this analysis.

## Required sequence

1. Resolve the target identity, species, construct and binding site. Request site
   clarification only if available evidence/user information does not resolve it.
2. Retrieve experimental PDB structures. Verify target/chain mapping, distinguish
   biological states and constructs, inspect local completeness and alternate
   locations, and deduplicate chains/depositions without losing source provenance.
3. Align a suitable stable protein core and compare the same pocket region/grid.
   Report cavity intersection/union, local subpocket changes and protein-derived
   donor/acceptor, hydrophobic, aromatic and charge features with typing provenance.
   Missing atoms must not create artificial pocket openings. Water/metal effects
   remain explicit evidence or unsupported features, never silently inferred.
4. Cluster pocket states using geometry and relevant chemical differences. Keep
   shared pocket cores; distinct states need not have disjoint volumes. Thresholds
   require validation, not invented universal cutoffs. PDB cluster membership
   denotes evidence support, not equilibrium population or occupancy probability.
5. Select real experimental medoids subject to structure quality. Retain supported
   rare states for review; small cluster size alone does not justify exclusion.
6. Show representative structures, overlap matrix, difference views, provenance,
   uncertainty and inclusion/exclusion reasons. Record explicit user adoption or
   edits, tied to a hash of the analyzed structure set and selected representatives.
7. Compute interaction recurrence within each adopted pocket state, distinguishing
   cross-state shared contacts from alternative state-specific contact patterns.
   Avoid counting duplicate structures as independent evidence. Never require a
   ligand to satisfy incompatible contacts from multiple receptor states at once.
8. Continue with state-aware 3D search, receptor-matched docking, N-E rescoring and
   later affinity models. Preserve pocket_state_id, receptor ID/hash and pose IDs.

One usable state can result in one representative after review; it is not evidence
that all possible states were observed. No usable PDB means this branch is blocked
with an explicit evidence gap. Cofolding and other generated states are deferred.

## Implementation and acceptance boundary

Planner text describes this requirement and accurately exposes the missing
capability. It is not a runtime gate and does not retrofit old sealed plans.
The future executor must require an owned, current adopted pocket-state receipt
before dispatching the new project's consensus; changing inputs invalidates it.

Before enabling that path, validate rigid-transform invariance, identical-pocket
agreement, local occlusion and feature-direction sensitivity, missing-atom handling,
duplicate suppression and per-state contact provenance. Compare representative
coverage with additional same-cluster structures on a fixed ligand panel. MDM2 is
the first intended real-structure acceptance case, not yet an executed result.

Method background: https://pmc.ncbi.nlm.nih.gov/articles/PMC5751414/
MDM2 state evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC4104591/
