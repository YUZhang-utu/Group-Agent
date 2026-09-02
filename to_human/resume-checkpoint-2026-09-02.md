# AIDD Macrocycle Agent — Resume Checkpoint

Date: 2026-09-02  
Repository: `https://github.com/YUZhang-utu/Group-Agent.git`  
Implementation baseline: `36b0f58`  
Validated local suite: `74 passed`

## Hard project boundary

Continue only `aidd_macrocycle_agent`, the Group-Agent AIDD macrocycle project.
Do not read, modify, or mix in the separate KRAS necessity/enhancement project.

Desktop repository:

```text
D:\agent\projects\aidd_macrocycle_agent
```

Linux workstation repository:

```text
/mnt/medchem_taltio/wrk/yu_agent/Group-Agent
```

## Current scientific objective

Build a fast, appendable 3D similarity workflow for eventual hundred-million to
billion-scale conformer libraries. Library-side chemistry is computed once;
new molecules enter as immutable artifact shards and receive stable int64 IDs.
Anchor evidence may rerank but must never reduce recall.

## Validated library and query

- Library `LIB-AFA68EE6888C`: 100,000 molecules, 299,999 conformers.
- Two real artifact shards: 149,999 and 150,000 conformers.
- Query: WEE1 PDB 8BJU, ligand QT9, chain A residue 601.
- 8BJU resolution 1.53 Å; ligand RSCC 0.933 and EDIA 0.951.
- Query chemistry uses authoritative QT9 CCD connectivity.

## Completed 3D infrastructure

- Int16 0.01 Å heavy coordinates and pharmacophore centers.
- Feature/meta/PMI/bounds/USRCAT artifacts, mmap offset reads, checksums, and
  append-only shard catalog.
- Complete 299,999-conformer artifact build: about 395 seconds with 24 workers.
- Incremental FAISS IVF-PQ m20 index with stable IDs and direct map disabled.
- QT9 final-index recall at nprobe 256: top-100/top-1,000 = 1.0;
  top-10,000 = 0.9969 in a 100,000-candidate pool; warm query about 14 ms.
- Candidate schema stores named objective-specific poses/transforms, unweighted
  and anchored atom-centered/projected scores, and raw per-anchor cross/self
  overlaps mapped through stable anchor IDs in `query_manifest.json`.
- Anchor, projected-feature, and pharmacophore-count information is strictly
  reversible reranking evidence. L1 recall remains anchor-independent.

## Anchor and Reduce status

- Hydrogen-independent QT9 evidence identified three polar contacts at
  2.876, 2.983, and 2.838 Å. Relative SASA burial is 1.000, 1.000, and 0.990.
- Windows has no Reduce; it computes distance/SASA only.
- Linux `aidd-workstation` has Reduce configured through an explicit `-DB` HET
  dictionary path.
- The generated CCD-to-legacy dictionary was not recognized by the installed
  Reduce and must not be used for QT9 production validation.
- The official wwPDB legacy QT9 HET entry succeeded. The user confirmed
  `accepted: true` and 34 ligand hydrogens added. Protein hydrogenation and a
  GLN375 flip were also reported.
- Pocket-fragment `appear unbonded` messages are expected truncation warnings;
  missing QT9 connectivity remains fatal.
- Code for Reduce-derived D-H...A angles and observed protein-partner projection
  points is implemented and pushed, but the real Linux enrichment output has
  not yet been run/reported.

## Immediate next action on Linux

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull origin main
python -m pip install --no-deps -e .

python -c "
from pathlib import Path
from aidd_agent.anchor_extraction import extract_query_manifest
m = extract_query_manifest(
    Path('data/e019_query_8bju/8BJU.cif'),
    Path('data/e019_query_8bju/QT9.cif'),
    'QT9', '8BJU:QT9:A:601')
print(m.write(Path('data/e019_query_8bju')))
"

aidd-agent enrich-query-anchors \
  --query-manifest data/e019_query_8bju/query_manifest.json \
  --ccd data/e019_query_8bju/QT9.cif \
  --hydrogenated-pdb data/e019_query_8bju/reduce-QT9-v3/query-pocket.reduce.pdb \
  --output data/e019_query_8bju/query_manifest.reduce.json
```

Adjust `reduce-QT9-v3` only if the successful run used another directory.
Inspect each anchor's `angle_degrees`, `hydrogen_acceptor_distance_angstrom`,
`angle_status`, donor hydrogen, and `projected_points`.

## Next development after real angle validation

1. Validate all three QT9 directional contacts and projection coordinates.
2. Implement Gaussian shape plus atom-centered/projected color overlap under
   named objective-specific rigid transforms.
3. Preserve raw per-anchor cross/self overlap for arbitrary later reranking.
4. Run shape-only vs unweighted color vs anchor/projected ranking diagnostics.
5. Add protein exclusion-grid filtering and terminal-torsion refinement.
6. Acquire Platinum/PDBbind conformer coverage data and LIT-PCBA enrichment
   benchmarks before promoting a native billion-scale kernel.

## Storage and source-control boundary

`data/`, `6aa/`, structures, indices, Reduce outputs, and machine-local profiles
are ignored and remain local. Git contains only code, tests, environment
templates, protocol, and research records. No new MOL2 batch is currently
required; use a later independent shard for final incremental acceptance.
