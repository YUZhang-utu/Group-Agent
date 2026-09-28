# Macrocycle conformer blocks: first offline baseline

Implemented entry points: `aidd_agent.macrocycle_source_probe` for bounded MOL2
inspection and `aidd_agent.macrocycle_blocks` for a uniform conformer sample from
the paired precompute catalogs. These are CLI tools; chat dispatch is not wired yet.
Neither entry point changes screening, source structures or existing registry IDs.

## Source inspection

The supplied `6aa/split_0001.mol2` and `split_0002.mol2` contain explicit hydrogens,
original molecule/conformer names and atom IDs. Inspected atom substructures are
`UNK`; some bonds carry `BACKBONE|DICT|INTERRES`. Names such as
`c--A-W-Lnme-A-dLnme-dVNMe-c_conf1` preserve building-block tokens, but no verified
token-to-atom correspondence has been established.

The final prefix probe reads 1000 records from each file. Of 2000 conformers,
1750 are assigned and 250 require an explicit ring mapping because another
perceived ring shares multiple atoms with the macrocycle. Assigned records have
18 ring atoms, six amides and six graph-derived peptide units; they cover 584
source molecule names. With diagnostic capacity 1000, seven blocks result.
This prefix is not a representative random library sample. Unassigned means
unresolved extraction, not chemically invalid or rejected from screening.

## Run and inspect

From the repository root in the existing scientific Python environment:

```powershell
python -m aidd_agent.macrocycle_source_probe --mol2 D:\agent\projects\aidd_macrocycle_agent\6aa\split_0001.mol2 D:\agent\projects\aidd_macrocycle_agent\6aa\split_0002.mol2 --output data/macrocycle-prefix-new --per-file 1000 --capacity 1000
```

For a precomputed library on the workstation:

```bash
python -m aidd_agent.macrocycle_blocks \
  --batch /path/to/library-precompute \
  --output /path/to/new-macrocycle-pilot \
  --limit 100000 --capacity 10000 --seed 20260925
```

The batch must contain `artifacts/catalog.json`, `chemical/catalog.json` and
`registry.sqlite3`. Catalog identity and the chemical-to-geometry catalog hash
must agree. Repeat in new directories with capacities 20000 and 30000 and the
same seed for capacity comparisons. The pilot limit is at most 100000 conformers.

Outputs:

- `report.json`: counts, scope, reasons and output hashes.
- `descriptors.jsonl`: original identity/provenance, ordered ring indices,
  torsions, amide states and graph-derived peptide units.
- `memberships.jsonl`: conformer-to-block mapping, retaining molecule identity.
- `blocks.jsonl`: block counts and descriptor centroids/radii.
- `unassigned.jsonl`: unresolved records and reasons; never silently discarded.

The direct source probe uses explicitly local IDs, not production registry IDs.
Its heavy-atom source-ID array translates zero-based descriptor indices back to
original MOL2 atom IDs. The catalog adapter preserves registry IDs and uses
zero-based heavy-atom artifact indices. Optional `--ring-mapping` JSONL records
require `conformer_id`, matching source `content_sha256`, and an ordered
`ring_atoms` list in that heavy-atom index basis. The list must form a closed
bond path. Supply mappings only after graph verification.

## Method and limits

Hard groups use ring size, cyclic atom/charge/bond tokens and amide states.
Absolute amide torsions up to 30 degrees are cis, at least 150 are trans;
intermediate values remain a separate boundary state. These are exploratory
geometric bins, not inferred energy barriers. Within each group, recursive median
splits on the largest-variance sin/cos coordinate enforce a maximum capacity.
Blocks are not equal-size chemical clusters; small groups remain small.

Cyclic traversal removes dependence on the starting atom and traversal direction.
Symmetric backbones use a torsion tie-break that still needs boundary-stability
benchmarking. Eligible backbone single-bond count is not independent ring freedom.
Building-block identity, side-chain properties/orientation, named phi/psi basins,
principal-axis geometry, incremental insertion and scalable full-library fitting
remain to be implemented. Different molecular identities remain separate records,
although chemically different side chains can share a backbone-only block.

Existing molecule-level block-tail utilities cannot directly consume multiple
conformer blocks per molecule; that adapter and molecule-level evaluation are
pending. No block priority, top-N enrichment, speedup or recall is claimed here.

## E076: supplied CSV identities (2026-09-25)

The two supplied CSVs each contain 50000 rows with `SMILES` and `Name` columns.
The bounded audit joined all 668 names in the previous 2000-conformer source
prefix to exactly one CSV row across both CSV files. All 2000 strictly sanitized
MOL2 structures match the SMILES heavy graph and specified stereochemistry;
this includes the 250 records with unresolved main-ring selection in E075.
This is chemical identity evidence, not a resolution of their main-ring maps.

`python -m aidd_agent.macrocycle_identity_audit --mol2 FILE1.mol2 FILE2.mol2 --csv FILE1.csv FILE2.csv --output NEW_DIRECTORY --per-file 1000`

Outputs are `report.json` and `identities.jsonl`, including CSV row provenance,
canonical structures and one atom correspondence where verified. Atom maps are
not asserted unique under symmetry and do not assign building-block labels.
Actual local results: `data/e076-csv-audit/`. Source files were not modified.
