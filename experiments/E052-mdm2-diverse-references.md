# E052: MDM2 crystal census and diverse reference selection

Protocol recorded before the census/similarity run. Exploratory project workflow.
Preserve D:/agent/MDM2 original files and the running E050 workstation task.

Question: which experimentally resolved human MDM2 ligand-bound structures
provide quality-qualified, chemically diverse references for full-library 3D search?

Inventory every supplied PDB/CIF with hashes, actual chain and residue content;
coordinate-only aligned exports cannot supply original crystallographic metadata.
Verify human MDM2 identity from UniProt and target-associated RCSB entries.
Retrieve all experimental entries without a resolution filter for the census.
Report distinct PDB counts in <=1.5, (1.5,2.0], (2.0,2.5], (2.5,3.0], >3.0,
and unknown resolution bins. Keep small-molecule, peptide/protein partner,
additive-only and method distinctions explicit. These categories can overlap.

Download original coordinates for all X-ray structures <=3.0 Angstrom, plus
user-provided IDs. Verify target-associated chains, ligand occupancy/alternate
states, heavy-atom completeness against CCD and target contacts. Quality default
for reference proposals is <=2.5 Angstrom, complete ligand heavy atoms, occupancy
>=0.9, no unresolved alternate conformer, and >=3 target residues within 4.5 A.
This is a quality screen, not experimental density validation. Preserve failures.
Use the supplied 5C5A pocket as the site reference after sequence mapping; require
at least three shared contact residues before mixing reference hypotheses.

Compare nonpolymeric organic ligands using RDKit Morgan radius-2, 2048-bit,
chirality-aware Tanimoto. Deduplicate identical canonical isomeric SMILES; retain
all crystal occurrences. Do not infer chemical graph diversity from coordinate
RMSD. Peptide/protein partners remain separately reported and are not silently
discarded or forced into small-molecule fingerprints.

Propose a deterministic farthest-first panel, seeded by the best-resolution
qualified occurrence and stopping when every qualified unique ligand has >=0.6
similarity to a reference, or 8 references are reached. Both settings are explicit
exploratory display defaults. Report uncovered ligands at the cap and a full
pairwise similarity matrix. Final reference acceptance remains a human review.

Outputs: raw API snapshots, downloaded mmCIF, local inventory, resolution census,
instance quality/contact table, ligand similarity matrix, proposed reference set,
and readable report. Build reusable code and a chat action for this survey.
This stage does not launch another full-library computation or claim that the
single-query guided adapter already performs multi-reference aggregation.
