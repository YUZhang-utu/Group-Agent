# MDM2 contact and spatial mapping audit

Date: 2026-09-23. Result: exploratory coordinate audit complete; grouped scoring
is not implemented or adopted. Production extraction and screening are unchanged.

## Evidence and reproducibility

Audited all 48 prepared ligand instances from 45 PDBs in the existing local
model-final cohort. Original source hashes passed. All prepared ligand shape
coordinates matched the saved protein transforms exactly (maximum discrepancy
0.0 A); pairwise distance invariance was checked for every complex. No ligand
atoms were removed by the audit's quality filter. Receptor atoms with alternate
locations, occupancy below 0.9 or insertion codes were excluded and counted;
unobserved geometry is not a negative contact observation.

The independent spatial probes are the observed p53 W23, L26 and F19 side chains
in local 1YCR, mapped to P04637 by explicit database correspondence and exact
peptide sequence. MDM2 maps to human Q00987. Alignment to the 5C5A reference frame
used 51 local MDM2 C-alpha pairs: RMSD 0.7893 A, maximum residual 1.4541 A.
This is an exploratory reference-volume comparison, not a new receptor-state
admission or validation of induced fit.

Protocol: [E053-subpocket-audit-protocol.md](../experiments/E053-subpocket-audit-protocol.md).
Script: [audit_mdm2_subpockets.py](../scripts/audit_mdm2_subpockets.py).

```powershell
$env:PYTHONPATH='src;D:/agent/tmp/e052-deps'
D:/conda/python.exe scripts/audit_mdm2_subpockets.py `
  --survey D:/agent/MDM2/analysis/e053/model-final/report.json `
  --peptide-reference D:/agent/MDM2/analysis/e052/structures/1YCR.cif `
  --output data/e053-subpocket-audit-new-run
```

Output directory must be new. The first attempt stopped before computation
because the existing receptor mapper rejects short peptide sequences; the audit
now verifies the peptide through exact mmCIF database/sequence correspondence
without changing that receptor admission policy. Successful output:
`data/e053-subpocket-audit-20260923-v2/`.

## Contact representation loses substantial evidence

Counts below are complexes, not independent molecules or validated interactions.
The atom proxy measures proximity between ligand C/S/halogen atoms and receptor
side-chain C/S/halogen atoms; it is broader than the current pharmacophore
extractor and is not a validated hydrophobic or halogen-bond assignment.

| MDM2 residue | Existing hydrophobic-mode support | Atom proxy <=4.0 A | Atom proxy <=4.5 A |
|---|---:|---:|---:|
| LEU54 | 25/48 | 48/48 | 48/48 |
| LEU57 | 3/48 | 36/48 | 48/48 |
| PHE86 | 0/48 | 29/48 | 43/48 |
| ILE99 | 2/48 | 46/48 | 48/48 |
| ILE61 | 39/48 | 48/48 | 48/48 |
| MET62 | 0/48 | 42/48 | 46/48 |
| TYR67 | 4/48 | 26/48 | 42/48 |
| VAL93 | 20/48 | 46/48 | 48/48 |

These gaps precede the LLM's selection of 20 modes. Concrete raw examples:

- 3JZK:YIN:A:1 ligand BR2 to PHE86 CE2: 3.7230 A, ligand atom maps to Trp23 probe.
- 6Q9L:HTZ:A:201 ligand CL1 to PHE86 CZ: 3.9475 A, ligand atom maps to Trp23 probe.

Code inspection explains coverage limitations: only feature type 5 (Hydrophobe)
feeds hydrophobic hypotheses; aromatic features use a separate branch; unknown
families including LumpedHydrophobe are omitted from the prepared feature table.
The receptor atom list for MET contains only CB. A raw aromatic/halogen atom
contact is therefore not guaranteed to become a hydrophobic mode. These are
representation limitations, not evidence that the observed contact is absent.

GLY58 has no side chain, so its side-chain proxy count is zero by definition.
Any proposed GLY58 contribution needs a separate backbone-contact definition.

## Spatial probes support the user's concern, but not residue-only groups

Each ligand heavy atom is assigned to its nearest p53 side-chain atom set only
if within the cutoff. A nearest/second-nearest difference <=0.5 A is ambiguous,
and one ligand atom can never occupy two probe groups. These are spatial overlap
proxies, not calibrated cavity occupancy rules.

| Probe | 1.5 A | 2.0 A | 2.5 A |
|---|---:|---:|---:|
| A: Trp23 | 48/48 | 48/48 | 48/48 |
| B: Leu26 | 46/48 | 46/48 | 46/48 |
| C: Phe19 | 47/48 | 47/48 | 47/48 |
| All three in the same crystal pose | 45/48 | 45/48 | 45/48 |

Thus 100% Trp23 spatial overlap is observed in this admitted cohort under all
three specified probe cutoffs. It is not universal binding necessity or a claim
about excluded structures, other receptor states or the screening library.

Exceptions requiring inspection, not automatic rejection:

- 4JV7:1MN:A:201 and 4JV9:1MO:A:201 have no Leu26-probe overlap.
- 4ZFI:4NJ:A:201 has no Phe19-probe overlap.
- 4JV9:1MO:A:201 is one of the user's ten recommended templates. A strict
  three-pocket gate would conflict with that reference under this definition.

The user's residue-group OR at 4.5 A reports A/B/C contact proxies in 48/48
complexes each, whereas the spatial probes report 48/46/47. Residue membership
alone would mask the spatial exceptions.

At the 2 A spatial cutoff and 4 A heavy-atom contact cutoff, ILE99 contacts ligand
atoms assigned to A in 39 complexes and B in 28. LEU54 contacts A in 48 and B in
45. VAL93 spans A/B/C (31/27/39); ILE61 spans A/C (46/43); GLY58 spans A/C (46/33).
These are overlapping complex sets, not additive counts. ILE103 has one A contact
and no B contact under this particular definition. The proposed residue lists
must therefore be spatially qualified rather than treated as exclusive labels.

## Polar-contact audit

The following counts indicate a prepared ligand donor feature within 3.5 A of
the named protein atom. They do not establish a hydrogen bond without chemistry,
direction, protonation and competing-partner review.

| Protein atom | Complexes with observed atom | Donor-proximity positives |
|---|---:|---:|
| LEU54 backbone O | 48 | 11 |
| GLN72 backbone O | 43 | 4 |
| GLN72 side-chain OE1 | 42 | 1 |
| VAL93 backbone O | 48 | 4 |

VAL93 proximity occurs in 4HBM:0Y7:A:201 (3.2439 A), 4ERE:0R2:A:201
(3.2466 A), 4QO4:35S:A:201 (3.1424 A), and 4OCC:2TZ:A:201 (3.4045 A), despite
no VAL93 HBD mode in the full consensus. The existing hydrogen-bond extractor
retains only the nearest compatible protein partner per ligand donor/acceptor;
these additional pairs are candidates for review, not confirmed missing bonds.

LEU54's current consensus mode has 10/45 distinct-PDB support and remains
ineligible as a mandatory anchor. Neither this coordinate audit nor donor
proximity measures full-library NH retention.

## Artifacts and next decision

- `report.json`: alignment, per-complex overlap at each cutoff, source/code hashes.
- `atom_contacts.csv`: every quality-observed heavy-atom pair <=4.5 A, with atom
  names, canonical residue, distance and spatial assignment of the ligand atom.
- `residue_contacts.csv`: minimum raw and nonpolar-proxy distances per residue.
- `spatial_atoms.csv`, `reference_probes.csv`: aligned coordinates and distances.
- `polar_distance_checks.csv`, `residue_summary.csv`: compact review tables.

Three focused geometry/quality tests passed: shared-boundary ambiguity, joint
rigid-transform invariance and rejection of uncertain atoms. Runtime checks
validated all 48 coordinate frames and source fingerprints.

Before grouped scoring: define chemically typed atom contacts independently of
the sparse feature table; review the three spatial exceptions and overlapping
residue memberships; retain separate alternative receptor states. Then implement
optional family maxima and nonlinear subpocket preference on one pose, preserving
underlying contacts. Do not introduce a hard three-pocket requirement or LEU54
mandatory rule from these results. Full-library NH and independent-active
retention remain separate pending checks.
