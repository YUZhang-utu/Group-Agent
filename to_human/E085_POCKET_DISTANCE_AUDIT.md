# E085 pocket distance audit

The original E084 distance matrix is superseded for scientific interpretation.
Four chemical channels were geometrically suppressed by a code defect. The
original six clusters are not an established set of MDM2 pocket states.
All E085 classifications remain exploratory and unadopted.

## Identity and exclusion

22IZ is authentic: the Q00987 RCSB entry search cache includes it, entries.json
describes it, and the CIF contains `_entry.id 22IZ`. The official record is
https://www.rcsb.org/structure/22IZ: MDM2 with St3-R stapled peptide, X-ray 1.52 A,
released 2026-08-26. This is not a correction to 2ZIZ. Exact source hashes and
the original query are exported in 01_22IZ_provenance.json.

7BJ0 is excluded before chain processing, at entire-entry scope, by explicit
user instruction. Original files are preserved. E085 contains 67 chains from
59 PDB entries. Every remaining cavity mask and grid coordinate is identical
to its E084 counterpart, verified by array equality. This isolates descriptor
changes from alignment changes. The original 7BJ0 alignment failure remains
a separate geometry finding, not an explanation of the inactive fields.

## Field defect and repair

The old nondirectional fields had radius 2 A from protein atom centres. The
cavity excluded C up to 2.20 A, N to 2.05 A and O to 2.02 A. Thus their fields
were buried entirely inside excluded space. Donor/acceptor points were shifted
outward and escaped this defect. This was a spatial-support bug, not absent
chemical labels or evidence that MDM2 lacks hydrophobic pockets.

Descriptor v2 defines nondirectional contact support from the atomic exclusion
surface outward by feature_radius (default 2 A); directional proxies retain
their prior definition. It is a residue-template contact annotation, not an
interaction energy or electrostatic calculation. Shell-width sensitivity holds
donor/acceptor fields fixed to isolate the changed channels.

After repair, hydrophobic and aromatic fields are nonzero in all 67 chains.
Positive is nonzero in 36, negative in 6, donor in 62, acceptor in 67. Remaining
zeros are explicitly checked against the nearest accessible grid distance:
no free grid point falls within that channel's support. They are local-region
zeros, not evidence that the whole protein lacks charged residues. See the
per-atom-source counts and support margins in tables 03 and 10.

The p53 Phe19/Trp23/Leu26 recognition subsites motivate these annotations;
these are p53 ligand residue labels, not MDM2 residue numbers. The seed region
is still reference-ligand expanded, not exhaustive access-path or cavity mapping.

## C1 reassessment and requested split

The old 55-member C1 is exported as two diagnostic strata: 43 legacy-zero and
12 legacy-nonzero members relative to 6Y4Q:A. Only one member, the representative
itself, has zero *uncensored* local difference. Therefore 42 of the 43 zeros
were produced by the hard support trigger, not identical local geometry.

We retained three local matrices:

- Legacy: maximum regional IoU loss only when changed volume >=30 A^3.
- Raw: maximum regional IoU loss without the trigger. Tiny accessible unions
  can have loss 1, so this strongly overfragments the collection.
- Smooth: maximum of regional loss * changed_volume/(changed_volume+30 A^3).
  This retains nonzero subthreshold differences, with no hard discontinuity.
  It is an exploratory regularization, not a validated universal metric.

For the proposed 43/12 split, average smooth within-stratum distances are
0.3195 and 0.3265; between strata is 0.3371. Corrected composite averages are
0.2397 and 0.2577 within, versus 0.2620 between. These overlap substantially;
the old zero/nonzero split does not establish two well-separated states.

The requested local-only reclustering is fully exported, with weighted average
linkage and unchanged one-unit-per-PDB weighting. At local cutoffs
0.20/0.25/0.30/0.35/0.40/0.45 it gives 26/16/8/4/2/1 candidate groups.
The two-group cut at 0.40 is **54+1**, not the old apparent 43+12 modes, and
its large group's representative radius is 0.416 >0.40. Do not adopt it as a
validated two-state partition. Table 04d contains every candidate membership,
representative and radius warning, including the four-group cut at 0.35.

Additional fixed complete-linkage diagnostic cuts (raw 0.50 and smooth 0.25)
are exported to demonstrate the fragmentation cost of bounding every pair.
Their many small groups are not proposed receptor families.

## Sensitivity and C6

All rows below exclude 7BJ0. These variants are not interchangeable calibrations.

| Descriptor | Cut 0.30 | Cut 0.40 | Cut 0.50 |
|---|---:|---:|---:|
| Old fields and old local trigger | 11 | 5 | 2 |
| Corrected fields, old local trigger | 9 | 4 | 2 |
| Corrected fields, uncensored local maximum | 56 | 31 | 10 |
| Corrected fields, smooth local support | 3 | 1 | 1 |

This separates the chemical fix from the local-metric change. It does not
validate a new one-cluster answer. At composite cutoff 0.30, the smooth support
scale 10/30/50 A^3 gives 13/3/2 groups. Shell width 1/2/3 A gives 3/3/2 groups
with directional fields held fixed. This is substantial parameter dependence.

C6 (8GCG:A) is retained in every calculation. With repaired chemistry and the
legacy local metric, it is singleton at 0.30 and joins 55 old C1 members at
0.40 (56 total). Under the smooth variant it belongs to a 58-member group at
0.30. Its independent-state status is not stable, so neither deletion nor
adoption as a distinct receptor state is justified by this audit.

Average-linkage cutoff bounds the merging average, not every member's distance
to the representative. Tables explicitly count and list radius exceedances.
A representative cannot be claimed to cover every member at the linkage cutoff.
If that coverage guarantee is required, it must be a separate clustering or
representative-selection constraint, with its fragmentation assessed explicitly.

Sensitivity on MDM2 does not remove MDM2 tuning bias. Target-independent
thresholds still require independent targets, alternative reference seeds,
grid-resolution checks and ligand accommodation/cross-docking validation.
None of these outcomes claims PLANTS, affinity or N-E validation.

## Deliverables and reproduction

Local package: `to_human/E085_POCKET_REVIEW_FINAL.zip`. The unpacked directory
contains provenance, exclusions, channel volumes/support margins, all 2211
pair comparisons, old-C1 subdivisions, 72 sensitivity settings and their full
memberships, distance matrices, scenario NPZs and corrected grids/report.
SHA256SUMS.json covers the package contents. Absolute source paths in reports
refer to this workstation; original CIFs remain in the previous E084 export.

```powershell
$env:PYTHONPATH='src'
python -m aidd_agent.pocket_states build --source PATH_TO_E052/report.json --output NEW_RUN --reference-query 5C5A:NUT:A:201 --target-chain A --exclude-entry 7BJ0
python scripts/audit_e085_pockets.py --old OLD_E084/report.json --new NEW_RUN/report.json --output NEW_AUDIT
```

The audit script intentionally checks this fixed MDM2 experiment's reference,
22IZ source evidence and unchanged masks. It is not a general-target CLI.
Old v1 reports are blocked from new adoption; recomputation is required.

Validation: 628 passed, 3 skipped before the final adoption guard fixture;
the final focused regression and artifact verification are recorded in the
E085 validation receipt. No user adoption or downstream consensus was executed.
