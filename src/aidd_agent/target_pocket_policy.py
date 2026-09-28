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
Implementation status: pocket-space/chemical-feature clustering and its persisted
adoption gate are not yet implemented. Existing structure_diversity and
structure_consensus do not satisfy this prerequisite. Explain this limitation for
new end-to-end project requests; offer supported identity/PDB evidence collection
without claiming the full sequence can execute. Do not invent a pocket action,
cluster assignment, approved representative or receipt. Existing explicitly
requested legacy analyses remain available but are not pocket-state-reviewed.
"""
