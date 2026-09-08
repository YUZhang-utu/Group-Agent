# Resume checkpoint: E028/E029 pre-docking refinement

Date: 2026-09-08

## Resume boundary

Continue the AIDD-agent pre-docking retrieval project from Git commit
`b401393`. Do not mix the separate AlphaFold 3 ligand-complex exploration or
the KRAS necessity/enhancement project into this workflow.

## Completed and accepted

- The real library contains 100,000 molecules and 299,999 conformers in two
  immutable E019 artifact shards.
- Incremental FAISS state, stable global IDs, offset-addressed coordinate and
  feature records, and restart-safe staged Gaussian scoring are implemented.
- QT9 FAISS validation recovered exact Top-1,000 fully and Top-10,000 at
  99.69%. Anchors only add provenance/reranking evidence and never remove the
  baseline recall lane.
- Real staged Gaussian timing on the workstation was 40.83 seconds for PCA
  coarse scoring of 299,999 conformers and 77.99 seconds for pair refinement
  of 7,704 conformers. Chunk hashes, deterministic IDs, and resume checks
  passed.
- E025/E026 multi-cocrystal aggregation is accepted for 8BJU/QT9 and
  1X8B/824: 12,834 union molecules, 955 shared molecules, 5,000 admitted
  molecules, and 5,834 receptor-specific docking tasks.
- E027 confirmed strong complementarity: zero shared Top-100 hits and only
  7.44% shared molecules in the union. Early retrieval must use a protected
  per-query union rather than intersection-only selection.
- E028 aligned point-cloud and pocket-exclusion QC is implemented and validated
  offline. Its geometry-only output is annotation, not a retrieval filter.
- E029 one-time chemical companion is implemented and validated offline. It
  stores atomic identity/charge/aromaticity/chirality, bonds/stereo, feature
  memberships, directional features, and bounded terminal torsions keyed by
  artifact-v1 global IDs.
- Directional Gaussian, element-aware soft vdW exclusion, chemical SDF export,
  and two-terminal-torsion reference primitives exist. The final batch-ranking
  integration and real workstation timing remain pending.
- Dependency-light test suite: 115 passed.
- Commits through `b401393` were pushed to `origin/main`.

## Raw-source boundary

`data/e019_artifacts/split_0001` and `split_0002` contain processed binary
artifacts only (`coords.bin`, `feats.bin`, `meta.bin`, `usrcat.f32.bin`, and ID
arrays). They do not contain the original MOL2 or enough information to
reconstruct authoritative atom/bond topology for E029.

The original files were confirmed on the Windows workspace:

```text
D:\agent\projects\aidd_macrocycle_agent\6aa\split_0001.mol2
D:\agent\projects\aidd_macrocycle_agent\6aa\split_0002.mol2
```

Expected SHA-256 values:

```text
bc729d337ba7065ebeccdbd7aa05f2985d4ecc778dd8b7949e1b57ada49a3b2b  split_0001.mol2
6458081fa50eb2676d818497403084775fa13fbc87cbcfbe6c7318f40141d1e7  split_0002.mol2
```

The user remembers uploading them to university storage, but their current
university path has not yet been located. If they cannot be found, transfer the
two files to this preferred workstation location:

```text
/mnt/local/hand/yuzhang/aidd/source/6aa/
```

Do not add these multi-gigabyte MOL2 files to GitHub.

## Immediate workstation continuation

1. Locate or transfer both original MOL2 files.
2. Verify both SHA-256 values before any real E029 build.
3. Pull `origin/main`, activate `aidd-workstation`, and install the editable
   package.
4. Run E029 with explicit source relocation overrides:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent

bash scripts/build_e029_chemical_companion.sh \
  /mnt/medchem_taltio/wrk/yu_agent/Group-Agent \
  --source split_0001=/mnt/local/hand/yuzhang/aidd/source/6aa/split_0001.mol2 \
  --source split_0002=/mnt/local/hand/yuzhang/aidd/source/6aa/split_0002.mol2
```

The expected companion catalog is:

```text
/mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json
```

5. Record wall time, worker utilization, peak RSS, shard reuse, conformer
   counts, output hashes, and any `.partial` state. This is the first real
   confirmatory E029 build; do not claim production throughput before it passes.
6. Run E028 with `CHEMICAL_COMPANION` set to the accepted catalog and inspect
   chemical SDF, vdW annotations, pocket contacts, and query/receptor frame
   identity.
7. Only after acceptance, connect directional and bounded terminal-torsion
   scoring into the fixed-budget final refinement stage and measure enrichment
   and latency deltas.

## Scale estimate and validation boundary

For 100 million molecules at three supplied conformers per molecule, the
existing 759-conformer/s artifact writer linearly extrapolates to about 4.6
days of core artifact construction. A realistic complete already-3D build
budget is approximately 1--2 weeks after indexing, verification, E029, and
recovery overhead. This is an estimate, not a measured 100M claim.

The intended warm per-query path remains fixed-budget FAISS/pharmacophore
union, capped PCA coarse scoring, and small pair/directional/torsion refinement.
Its provisional target is 2--10 minutes per co-crystal query. Before reporting
that target, physically validate 10M and 100M metadata/artifact stages, cold and
warm cache latency, shard-local Top-K recall, and interrupted resume.

If three-dimensional conformers must be generated from 100 million SMILES,
conformer generation becomes a separate weeks-to-months distributed-compute
problem and is not covered by the artifact-ingestion estimate.

