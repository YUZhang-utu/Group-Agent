# KRAS G12D/G12V: shared pocket contacts and allele-specific constraints

2026-10-23 exploratory coordinate inspection. This is a reference comparison,
not a finalized screening policy or evidence of universal interaction necessity.
The user requested agreement on necessary contacts before full-library submission.
No real KRAS library scan has run. The later user-confirmed shared-region
implementation and local validation are described in E116_KRAS_WEEKEND_RUN.md.
The comparison below records the evidence behind that decision.

Sources: https://www.rcsb.org/structure/7RPZ and
https://www.rcsb.org/structure/12FF. Coordinates downloaded from files.rcsb.org.
Local raw measurements: data/kras-weekend-20261023/contact_comparison.json.
Compare GDP-state native ligand/receptor pairs, chain A in each structure.

| Region | G12D 7RPZ / MRTX1133 | G12V 12FF / compound 20 | Interpretation |
|---|---|---|---|
| His95 | NE2 to ligand N01: 2.834 A | NE2 to ligand N4: 2.988 A | Similar N-acceptor/polar geometry in both references |
| Asp69 | OD1 to ligand O02: 2.630 A | OD1 to ligand N5: 2.642 A | Shared polar anchoring region, different ligand functional groups |
| Glu62 | OE2 to ligand N06: 2.847 A | OE2 to ligand N3: 2.701 A | Shared acidic contact region |
| Arg68 | NE to ligand N03: 3.010 A | Closest atom pair CB to ligand N5: 3.407 A | Same region, not the same polar contact mode; do not call the G12V pair an H bond |
| Met72 | CB to ligand C21: 3.598 A | CG to ligand C25: 3.646 A | Similar hydrophobic environment |
| Residue 12 | Asp12 OD1 to ligand N05: 2.703 A | Closest pair is Val12 backbone N to ligand C3: 4.262 A | A G12D acidic-sidechain polar constraint cannot be copied to Val12 |

Distances are heavy-atom coordinate measurements, not angle-validated hydrogen
bonds, energetic contributions or proofs of necessity. Both references contact
Tyr64/Tyr96/Gln99 regions as well, with different atom pairs and contact modes.
Arg68 is not established as a compulsory direct H bond for every G12V ligand.
12FF does not support requiring an Asp12-like polar interaction at Val12.

Sequence audit uses KRAS4B P01116-2. Both experimental constructs include C51S,
C80L and C118S backgrounds in addition to the requested codon-12 mutation.
7RPZ has an additional observed N-terminal numbering/identity discrepancy to
inspect. These are native crystallographic constructs, not silently repaired
canonical mutant structures. Their contact patterns must not be presented as a
general mutation-selectivity comparison across unrelated ligand chemotypes.

Implementation issue discovered during audit: RDKit BaseFeatures marks some
ligand fluorines as Acceptor features. Whole-color replay must preserve its
published feature definition, but an organic fluorine feature must not satisfy
a mandatory conventional N/O hydrogen-bond condition merely because its feature
type equals Acceptor. The delivered implementation now applies explicit N/O
member-atom eligibility and direct candidate/receptor proximity; this is not
angle-validated hydrogen-bond detection.

Recommended design for discussion: shared candidate decoding and cheap invariant
requirements, followed by independent native-state pose searches and acceptance.
G12D-specific branch can require the user-requested Asp12/His95/Arg68 pattern;
G12V branch needs Val12-compatible sterics and independently calibrated contacts.
Shared prefilter logic must be the UNION of branch eligibility. A candidate that
could pass either branch must survive. Never require a molecule to pass both
targets unless explicitly seeking dual-target hits. Shared I/O does not imply
twofold speedup: different reference poses and geometry still require computation.

Before submission: finalize interaction eligibility; test same-pose constraints,
negative controls, chunk integrity/resume and export; obtain workstation timing
and pass rates. Existing self-reference recoveries alone are insufficient.
