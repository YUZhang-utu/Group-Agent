"""Planning requirements, not an implemented pocket-clustering adapter."""

TARGET_POCKET_POLICY = """
For a new target-based screening project, require PDB pocket-state review before
interaction consensus or downstream screening. The intended order is verified
target/site identity, experimental PDB collection and quality control, aligned
pocket-space comparison with chemical features, receptor selection advice and real
structure representatives, explicit user adoption, then local interaction consensus.
Use experimental PDB structures in this first version; do not automatically add
cofolding predictions, MD frames or generated receptor conformations.
Ligand chemical diversity and global protein RMSD are not pocket-state clusters.
PDB counts are structural support, not equilibrium occupancy probabilities.
Do not automatically discard small clusters or equate high global pocket overlap
with unchanged local subpockets or interactions. Preserve state-specific contacts;
do not turn mutually exclusive state contacts into simultaneous mandatory anchors.
Implementation: the trusted coordinator runs pocket_states from structure_diversity,
then pocket_adopt on explicit user acceptance, then pocket_consensus for one adopted
pocket_state_id. Use pockets with reference {} or actual reference_query/target_chain
to prepare analysis. receptor_advice is the default user-facing result: provisional
single receptor for limited evidence or unresolved variation, stable-partition
representatives when supported, otherwise coverage representatives and recommended k.
Scaffold-series counts are proxies, not verified independent chemical series.
Use cited target-specific literature via literature_search for additional ideas;
distinguish published evidence, computed observations and hypotheses. Never infer
rigidity from few structures or physical continuity from an unstable partition.
The user chooses the receptor set. Show recommendation and reasons, not raw tables.
cluster_distance remains an exploratory diagnostic parameter;
Use adopt with optional pocket_selection {cluster_ids: actual reviewed IDs} to
adopt a subset or individual receptor_options; omit pocket_selection to explicitly
adopt the recommended set. Coverage neighborhoods are not discrete physical states.
inspect sensitivity counts instead of forcing a cluster count. Chemical fields use
residue templates and approximate directions, not full electrostatics or hydrogen
bond energetics. Never invent cluster IDs or adoption receipts. Model-generated
plans cannot emit source_run actions; these remain trusted coordinator operations.
New diversity outputs cannot bypass adoption via legacy structure_consensus.
Existing explicitly requested old analyses remain available without claiming review.
"""
