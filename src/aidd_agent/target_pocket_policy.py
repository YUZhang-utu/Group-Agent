"""Planning requirements, not an implemented pocket-clustering adapter."""

TARGET_POCKET_POLICY = """
For a new target-based screening project, require PDB pocket-state review before
interaction consensus or downstream screening. The intended order is verified
target/site identity, experimental PDB collection and quality control, aligned
pocket-space comparison with chemical features, pocket-state clustering and real
structure representatives, explicit user adoption, then state-specific consensus.
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
to prepare analysis. cluster_distance defaults to a permissive exploratory .40;
Use adopt with optional pocket_selection {cluster_ids: actual reviewed IDs} to
adopt a subset; omit pocket_selection to adopt all reported states explicitly.
inspect sensitivity counts instead of forcing a cluster count. Chemical fields use
residue templates and approximate directions, not full electrostatics or hydrogen
bond energetics. Never invent cluster IDs or adoption receipts. Model-generated
plans cannot emit source_run actions; these remain trusted coordinator operations.
New diversity outputs cannot bypass adoption via legacy structure_consensus.
Existing explicitly requested old analyses remain available without claiming review.
"""
