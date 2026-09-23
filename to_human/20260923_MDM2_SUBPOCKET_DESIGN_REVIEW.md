# MDM2 subpocket design review - 2026-09-23

Status: user-proposed design under review, not adopted or implemented. The user
has now reported a completed live consensus recommendation. No library execution
or adoption of the following edits is confirmed.

## Requested design

- Conditional mandatory contact: pocket:db6f8400e16ed844, ligand HBD to LEU54
  backbone O. First measure library NH retention; do not enable the rule yet.
- Subpocket A (Trp23): hydrophobic contacts to LEU57/PHE86/ILE99.
- Subpocket B (Leu26): LEU54 side chain/VAL93/ILE99/ILE103.
- Subpocket C (Phe19): ILE61/MET62/TYR67/VAL75/GLN72/GLY58.
- Occupancy: at least one validated member contact in a subpocket, evaluated in
  a single pose. Three occupied subpockets should outrank two strongly, then one.
  The user has not specified numerical nonlinear coefficients or a hard minimum.
- Optional grouped bonuses: HIS96 pi alternatives 0.60; HIS96 HBD/HBA
  alternatives 0.55; GLN72 backbone O hydrogen bond 0.50; VAL93 backbone O
  hydrogen bond 0.50; TYR67 contact 0.35 with a proposed four-finger annotation.
- Four previously recommended water bridges: weight at most 0.20 or removed;
  exact choice pending. Merge ILE61 hydrophobic modes into one scoring family
  with weight 0.40, preserving each underlying spatial mode.

## Inspection findings

Local evidence: D:/agent/MDM2/analysis/e053/model-final/report.json, not a newly
retrieved workstation report. The inspected IDs match the user's recommendation.

- consensus_model.consensus retains low-frequency modes. Its >=0.5 frequency,
  >=3 PDB and >=2 chemotype rule controls mandatory eligibility, not retention.
  Recommendation subsequently selects at most 20 modes.
- LEU54 pocket:db6f8400e16ed844 has support 10/45 and is not mandatory-eligible.
  Human edits are subject to the same gate; making it mandatory needs an explicit
  reviewed policy change as well as NH/geometry/positive-retention checks.
- Existing hydrophobic modes include LEU57 (one mode, three PDBs) and ILE99
  (two modes, two PDBs), but no PHE86 mode. The LEU57/ILE99 observation union
  comprises five query instances. This is NOT a validated occupancy frequency.
- GLN72 already has four backbone-O HBD modes plus one side-chain-OE1 mode:
  f079a27a67e5a4e9, bfa03273a49a2c5f, ddefa1b1c90b6290,
  138f4735025c6571 (backbone), 4af06df0fc296afe (side chain).
  VAL93 has no HBA/HBD mode in this report.
- Hydrophobic extraction uses ligand feature type 5, a 4 A distance, and nearest
  atom per residue. It is a feature-based geometric hypothesis extractor, not an
  exhaustive ligand-heavy-atom contact map. Missing modes need extraction and
  geometry review before concluding that a subpocket is unoccupied.
- Current alternative_groups are mandatory OR groups; they do not implement
  optional grouped bonuses. Ranking currently has weighted individual optional
  scores, not nonlinear subpocket occupancy scoring.

## Required semantics and verification before adoption

Retain atom/direction-specific evidence underneath functional subpocket groups.
Compute group support by a deduplicated union per evaluable complex, with missing
geometry reported as unknown rather than absent, and report PDB/chemotype support
separately. Do not sum individual frequencies or assume 100 percent occupancy.
Use aligned spatial regions as well as residue identities: ILE99 appears in A and
B, so a single ambiguous hit must not imply occupation of two distinct regions.
Do not pool contacts from different poses, templates or receptor states.

Optional mode families should count once (for example, maximum valid member
score), rather than earning a sum merely because several equivalent modes were
extracted. Preserve all member modes; do not average their coordinates.
Determine whether strong three-pocket preference means a nonlinear scalar score
or lexicographic occupancy priority before choosing coefficients. Avoid double
counting hydrophobic family bonuses and the occupancy reward without review.

Measure library NH presence by unique molecule and scaffold, not conformer count;
check original chemistry against prepared features to distinguish absent NH from
lost hydrogen/donor annotation. NH presence alone does not prove pose feasibility
for LEU54. Independent active retention is a separate check. No authoritative
full-library NH statistic was located in the inspected local records.

GLN72 O versus OE1 must remain distinct. A TYR67 contact in one fixed structure
does not establish induced fit or a four-finger binding mode. These annotations
require explicit structural evidence; new receptor states are outside the current
single-admitted-state path.

Primary structural reference located during review:
https://pmc.ncbi.nlm.nih.gov/articles/PMC11095393/
It describes LEU54 backbone hydrogen bonding and Trp23-site contacts in a peptoid
complex; it does not establish the proposed cohort's group occupancy frequency.
