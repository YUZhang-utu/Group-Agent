# Full source audit before macrocycle blocking

The audit reads every row and conformer in the explicitly supplied inputs. It is
not the previous bounded prefix probe. It does not infer coverage of additional
workstation shards from the two local shards.

## Completed local evidence, 2026-09-25

Full run: `data/e077-full-audit-v2/report.json`, 420.8 seconds with four workers.
Both CSVs contain 50000 rows; there are 100000 distinct names and no duplicate
canonical stereochemical structures among those names. All CSV residue sequences
validate, including their D/L and backbone N-methylation labels.
MOL2 counts are 149999 and 150000. All 299999 conformers pass chemical identity,
stereochemistry and peptide-unit mapping, with an 18-atom main ring. No CSV name
lacks a MOL2 record and no conformer name is duplicated. E075's previously
unresolved proline-containing cases are handled by peptide connectivity here.

Findings requiring explicit representation:

- `c--A-Wnme-PdFnme-dW-dLnme-c` has conf1 and conf2 only; conf3 is absent from
  the supplied sources. Do not invent the missing conformer or discard the molecule.
- 27658 names have different cis/trans/boundary patterns among their conformers.
  Excluding ambiguous cyclic alignments and boundary-only differences, 16653
  names show a definite cis-to-trans difference at the same named peptide bond.
- 16171 conformers have at least one boundary omega; 21212 have at least one cis
  omega. These sets overlap. Counts use the documented exploratory angular bins.
- 18 conformers have multiple equivalent name-to-cycle rotations, such as
  `c--A-Wnme-Fnme-A-Wnme-Fnme-c`. Preserve allowed rotations rather than assigning
  chemically meaningful absolute residue numbers to an arbitrary starting atom.

Example: `c--A-W-Lnme-A-dWnmedFnme-c` has omega at the fifth named peptide bond
of approximately 165.2, -14.3 and -158.8 degrees in conf1, conf2 and conf3,
respectively. Whole-molecule stereochemical identity matches in all three.
This directly supports conformer-level omega grouping while retaining one
original molecule identity. Additional query results are in
`data/e077-full-audit-v2/additional_checks.json`.

Validation: 28 focused tests passed, including explicit-H fragment handling,
proline cycles, D/L and methylation mismatches, symmetry and missing-source
accounting. The interrupted first attempt is explicitly superseded and excluded.

## Execution

Optional-CSV directory mode (E078): place any available CSV beside its MOL2 with
the same stem, for example `split_0001.mol2` and `split_0001.csv`. Other MOL2 files
may have no CSV. Directory traversal includes subdirectories and freezes the
discovered input list before execution. Keep the output outside the source tree.

```bash
python -m aidd_agent.macrocycle_full_audit \
  --source-dir /path/to/full-library-sources \
  --output /path/to/new-library-audit \
  --workers 4
```

Directory mode automatically allows missing CSV entries. Explicit file-list mode
remains strict unless `--allow-missing-csv` is specified; `--csv` is optional with
that flag. Available CSV rows are joined by exact molecule name across inputs,
not by row position. Unrelated differently named CSVs are not auto-discovered.

Each conformer records `verification_source` as `csv_and_mol2` or `mol2_only`,
plus the separate `csv_identity_verified` flag. MOL2-only acceptance validates
the supplied atom/bond graph, stereochemistry and coordinates; it does not infer
missing chemistry from coordinates alone or claim independent CSV verification.
Duplicate, invalid or conflicting CSV rows remain review cases, never a silent
fallback. The report also counts same-name/different-chemistry records.

If a MOL2-only name cannot be validated, graph-only peptide-cycle mapping can
retain geometry with `name_mapping=unverified`. Such rows have no established
cross-conformer residue alignment and are excluded from the corresponding
cross-conformer pattern-change statistic. Unsupported/non-peptide topology,
invalid chemistry or unspecified stereochemistry remain in the review queue.
Missing CSV by itself does not discard a conformer. Block admission must inspect
the audit flags, not assume that audit completion means all rows are validated.

E078 validation: 30 focused tests passed. A real mixed-source fixture of 24
conformers (12 with CSV, 12 without) passed; all 12 MOL2-only conformers retained
ring mapping and omega angles. This validates the new path, not the unseen full
workstation library. Evidence: `data/e078-mixed-audit/report.json`.

Run from the repository root in the scientific Python environment:

```powershell
python -m aidd_agent.macrocycle_full_audit --csv D:\agent\projects\aidd_macrocycle_agent\6aa\split_0001.csv D:\agent\projects\aidd_macrocycle_agent\6aa\split_0002.csv --mol2 D:\agent\projects\aidd_macrocycle_agent\6aa\split_0001.mol2 D:\agent\projects\aidd_macrocycle_agent\6aa\split_0002.mol2 --output data/new-full-audit --workers 4
```

Supply every intended CSV and MOL2 path for a larger-library audit. An omitted
shard is not audited. The SQLite name index and bounded worker queue avoid loading
the raw library into memory. Source hashes and per-record hashes are recorded.
Runtime and disk requirements still grow with library size; this is not a claimed
billion-molecule performance result. Interrupted runs require a new output
directory; resumable shard scheduling is not implemented in this entry point.

`report.json` updates during execution and says `complete` only after both scans
and cross-file accounting finish. `audit.sqlite` stores original names, CSV
provenance, molecular signatures, verified peptide units and conformer records.
`issues.jsonl` records failures and missing names; suffix counts and repeated
names are also reported. Completion means audit completion, not zero anomalies.
These source-level records do not replace production registry identities.

## Checks and interpretation

- Strip only a terminal conformer suffix for CSV name lookup. Never join by file
  row position; names may occur in different orders.
- Parse A/W/L/F/V/P/G with optional d and nme/NMe, including omitted hyphens.
  Unrecognized alphabets remain explicit review cases.
- Validate the entire ordered cyclic residue sequence against molecular graph
  fragments, including stereo and N-methylation. Reference fragments come from
  amino-acid structures, not a universal D=R or L=S lookup.
- Find the main ring by the directed peptide connectivity N-CA-C(=O)-N. Proline's
  own ring does not invalidate a uniquely defined cyclic peptide backbone.
- Cap cleaved endpoints only on analysis copies when comparing residues. This
  avoids different implicit-H conventions between SMILES and MOL2. No source
  hydrogen or stereochemical information is silently repaired.
- Check every conformer's whole chemical graph and stereochemistry, finite
  coordinates, original atom IDs, geometric omega state and name correspondence.
- Report missing molecules, duplicate names, duplicate conformer names and any
  deviation from the expected suffix set conf1/conf2/conf3. Such deviations are
  evidence for review, not permission to discard the affected molecules.

Omega uses CA(i)-C(i)-N(i+1)-CA(i+1). The stored state is cis for |omega| <=30,
trans for |omega| >=150 and boundary otherwise. Name prefixes d indicate residue
stereochemistry; they do not encode this amide geometry. A single molecule name
may therefore have different conformer omega patterns. Preserve identity and
record those patterns separately. Repeated identical sequences can produce more
than one valid cyclic alignment, which is explicitly counted.

## First blocking stage after audit

1. Keep a full manifest with one entry per original conformer and a separate
   unresolved queue. Join source identities to the actual registry before
   producing production block labels.
2. Use the verified peptide main-ring map, ring size and ordered omega pattern
   as the structural starting point. Preserve D/L, N-methylation, proline positions
   and residue identities for chemically valid cyclic alignment.
3. Compute phi/psi/omega sin/cos descriptors; record constrained units and eligible
   bond counts without claiming independent closed-ring degrees of freedom.
4. Compare capacities 10000/20000/30000 conformers within the same hard groups.
   Allow residual small blocks; do not mix incompatible states to fill a quota.
5. Add side-chain identity/properties and typed spatial features as separate
   benchmark variants. Do not silently promote phi-sign or psi bins to measured
   energy barriers.
6. Evaluate sampling 100/500 molecules per block using matched upper-tail fractions,
   with held-out molecule recall against the completed screening ranking. Reserve
   random exploration and deduplicate final outputs by molecule.

The full audit supplies validated mappings and exceptions. Full-library clustering,
incremental block insertion, chat dispatch and adaptive block exploration still
need implementation/integration and benchmarks before production use.
