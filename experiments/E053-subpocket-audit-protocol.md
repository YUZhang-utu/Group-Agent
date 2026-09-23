# E053 subpocket contact audit protocol

Date: 2026-09-23. Status: protocol written before execution.

Question: Are missing Trp23-site contacts absent from crystal coordinates, or
omitted by feature extraction/recommendation? Does residue-only grouping agree
with an independently defined spatial subpocket reference?

Scope: existing 48 prepared instances in model-final/report.json, with original
mmCIF structures, canonical human MDM2 mapping and saved protein transforms.
No new library scan, score modification, design adoption or model API call.

1. Verify original source hashes. Enumerate all quality-observed ligand and
   target-chain heavy atoms (model 1, no alternate location, occupancy >=0.9,
   no insertion code). Record excluded atoms; missing geometry is unknown.
2. Compute exhaustive ligand-heavy-atom to residue-heavy-atom distances and
   contact pairs <=4.5 A, independently of pharmacophore extraction. Compare
   distance thresholds 4.0 and 4.5 A; neither establishes interaction energy.
   Separately label carbon/sulfur/halogen proximity as a geometric proxy, not a
   chemically validated hydrophobic bond.
3. Verify saved coordinate transforms against prepared ligand shape points and
   check distance invariance. Preserve original evidence and explicit atom IDs.
4. Use local 1YCR as an independent spatial reference: verify its MDM2 Q00987
   chain and p53 P04637 sequence mapping; align the MDM2 local backbone to the
   5C5A reference frame using the existing cohort's local alignment residues.
   Report fit quality and side-chain atom coordinates for p53 W23/L26/F19.
   This defines exploratory spatial probes, not validated cavity volumes or
   admission of another receptor state into the screening cohort.
5. For each aligned ligand heavy atom, measure distance to each probe's nearest
   side-chain heavy atom. Assign to only the nearest site, leaving ties within
   0.5 A ambiguous. Report spatial-overlap hypotheses at 1.5, 2.0 and 2.5 A;
   missing occupancy under this probe definition is not absence of binding.
   Compare atom-based spatial membership with the user's residue groups and
   record which protein contacts belong to each spatially assigned ligand atom.
6. Compare residue-level raw contact support with existing full consensus
   observations and separately from the 20-anchor recommendation. Audit GLN72
   backbone O vs OE1, VAL93 backbone O and LEU54 backbone O without inventing
   hydrogen directions or assigning H bonds solely from heavy-atom distance.

Outputs: protocol, reproducible script, source/code hashes, atom contact CSV,
per-complex spatial assignments, residue summary and human review. Report both
complex and distinct-PDB counts. Any denominator must identify its evaluated
cohort; no extrapolation to all structures, chemotypes or library molecules.

Interpretation: exploratory audit. Confirmation is limited to deterministic
claims about the saved coordinates/extractor; occupancy thresholds and contact
chemistry remain hypotheses. Findings decide what extraction/mapping changes
are needed before grouped scoring can be implemented.
