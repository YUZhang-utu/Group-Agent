# E115 contact-region and feasible-pose review

Exploratory protocol, 2026-10-20. Implement three diagnostics without changing
production gates or launching a library search. Reuse the sealed E114 proposal
panel and verified E113 receipt through the E112 native-state preparation.

1. Review MOL-1EC0C7A6D657 (4382601) in its paired 6Q9L state. Export topology
   and coordinates, typed Gaussian overlap attribution, and isolated versus
   protein-occluded polar atom SASA. SASA is geometric exposure, not permeability
   or binding evidence. Feature membership defines polar atoms; report atom IDs.
2. Define groups using reference information only: identical sets of residue,
   protein atom and interaction-class labels, identical feature type and direction
   kind, complete-link distance <=2 A and directional agreement within 30 degrees
   (axial sign ignored). Each reference index belongs to exactly one group.
   Match groups one-to-one to candidate features using maximum member score.
   Report individual groups and original coverage in parallel. No new pass gate.
3. Start from both original and E114 refined endpoints. Retain all endpoints and
   all feasible trials; select the best feasible shape, preserving any original
   full-policy positive as a separate incumbent. Exact original contact, steric,
   envelope and region constraints determine feasibility. Deterministic bounded
   rigid coordinate search uses 0.20/0.08 A translation and 0.08/0.03 rad rotation,
   at most 96 trial evaluations per pose, and a cooperative total budget.
   A local result is not a global optimum; infeasible starts are not rescued.

Evaluate the 65 saved proposals and nine supported identity controls. Report
partial completion explicitly. Default budget 180 seconds excludes setup,
mandatory endpoint replay and primary-candidate export. Count poses, not unique
molecules or prevalence. Validate replay scores within 1e-6 and hash lineage.
Use synthetic unit tests for grouping, uniqueness, SASA and feasible retention;
workstation chemical evidence remains pending until the new receipt is returned.
