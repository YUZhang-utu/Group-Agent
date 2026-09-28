# Macrocycle audit and blocking handoff: 2026-09-26

## User decisions

- Partition actual 3D conformers, not molecules. conf1/conf2/conf3 may enter
  different blocks. Preserve original identities; deduplicate by molecule only
  for final docking delivery.
- CSV Name is the molecule name. MOL2 appends a terminal _confN suffix. Lowercase
  d denotes D amino acids; nme/NMe denotes backbone N-methylation. Verify names
  against structures, including names with missing hyphen separators.
- Audit the full library before production blocking. CSV coverage may be partial;
  available CSVs will be adjacent to their MOL2 with the same filename stem.
  Missing CSV must not by itself discard conformers. Use supplied MOL2 atom/bond
  graphs and coordinates and record the lack of independent CSV verification.
- Compare capacities 10000/20000/30000 conformers, allowing small residual blocks.
  Later benchmark backbone, chemistry and typed side-chain spatial variants.
  Adaptive exploration will sample 100/500 molecules per block and compare upper
  tails, with molecule-level recall against the completed screening ranking.

## Implementation and evidence

Repository: D:/agent/projects/aidd_structure_guided_chat.

E075: macrocycle_blocks.py implements an offline, bounded backbone torsion pilot.
It uses generic ring perception, geometric amide bins and median splits of sin/cos
descriptors. This baseline is not full production clustering. Initial prefix
ambiguities motivated the graph-verified peptide mapping below.

E076: macrocycle_identity_audit.py matches exact names against CSV and compares
chemical/stereochemical identity. Initial 2000 conformers / 668 names all matched.

E077: macrocycle_peptide.py verifies residue units and cyclic peptide connectivity,
including proline side rings. macrocycle_full_audit.py streams full supplied
sources using a disk-backed SQLite name index and bounded worker batches.
Completed full local audit in data/e077-full-audit-v2:

- 100000 unique CSV molecules, 299999 MOL2 conformers; all pass implemented graph,
  stereo and residue-name checks; all main rings have 18 atoms.
- 99999 molecules have three conformers. c--A-Wnme-PdFnme-dW-dLnme-c has only
  conf1/conf2; conf3 is absent from the two supplied files.
- No duplicate names, duplicate conformer names, missing CSV-to-MOL2 molecules
  or duplicate canonical stereochemical structure groups in these inputs.
- 16653 uniquely aligned names have definite cis/trans changes between their
  stored conformers. 16171 conformers contain boundary omega angles, and 18 have
  equivalent cyclic sequence rotations. Preserve these flags; do not rename or
  discard molecules to resolve them.
- Full input/output and implementation hashes are in report.json. Additional
  analyses are in additional_checks.json. Runtime was about 420.8 s, four workers.
- The first attempt data/e077-full-audit was interrupted after an analysis-fragment
  H-capping bug. Its superseded.json explains why its mismatch counts are invalid.
  Never use that first attempt as evidence of source chemistry errors.

E078: directory discovery and optional CSVs are implemented. --source-dir enables
recursive MOL2 discovery and same-stem CSV lookup; explicit file lists require
--allow-missing-csv for the fallback. Rows record csv_and_mol2 versus mol2_only and
csv_identity_verified. Existing CSV conflicts are review cases, never silent
fallbacks. Names that cannot be verified may retain graph-only peptide geometry,
with explicit alignment limitations. Unsupported chemistry remains a review case.
Mixed real fixture data/e078-mixed-audit has 24 passing conformers, 12 CSV-backed
and 12 MOL2-only. Latest focused suite: 30 passed; English guard and whitespace
check passed. Source data and production registry/search were not modified.

## Next action

1. Obtain the complete workstation source directory and confirm CSV placement.
   Only two local shards are available in D:/agent/projects/aidd_macrocycle_agent/6aa.
   No full workstation-library audit has been launched.
2. Review/deliver the current working-tree changes before asking the workstation
   to pull them. Changes from E075 onward have not been committed or pushed in
   these turns. Preserve unrelated build/ and to_human/E066_CHAT_UPDATE.patch.
3. Run the full audit with a fresh output directory outside the input tree:

```bash
python -m aidd_agent.macrocycle_full_audit \
  --source-dir /path/to/full-library-sources \
  --output /path/to/new-library-audit --workers 4
```

4. Check discovered manifest coverage, per-source hashes/counts, verification-source
   counts, missing/duplicate records and same-name chemistry conflicts. Missing
   CSV is acceptable; a completed audit does not imply zero anomalies.
5. Connect verified peptide maps to production registry identities and conformer
   descriptor extraction. The E075 generic-ring pilot has not yet been upgraded
   to consume the E077 peptide mapping throughout its production path.
6. Implement/benchmark full conformer blocking and molecule-aware exploration.

## Remaining work and limits

No production full-library clustering, scalable frozen cluster fitting, incremental
insertion, typed side-chain geometry, chat dispatch or adaptive exploration adapter
has been completed. Audit interruption currently requires a new output directory;
resumable shard scheduling is not implemented. Whole-library runtime and disk
requirements still scale with input size. Do not claim billion-molecule readiness.
Graph-only unnamed cycles lack verified cross-conformer residue alignment.
Omega bins (absolute <=30 cis, >=150 trans, otherwise boundary) are geometric
conventions, not measured energy barriers. phi-sign/psi-basin hard splits remain
benchmark hypotheses. Existing molecule-level tail utilities need adaptation for
conformers of one molecule occupying multiple blocks.

## Resume instructions

Read AGENTS.md, research-state.yaml, this handoff, the latest research_log.md entries,
to_human/E077_FULL_SOURCE_AUDIT.md and experiments/E078-optional-csv-audit.md.
Inspect git status before editing; preserve unrelated changes. Continue from
optional-CSV full-library audit preparation, not by rerunning the completed local
two-shard audit without a new reason.

Local verification environment:

```powershell
$env:PYTHONPATH='D:\agent\tmp\e052-deps;D:\agent\tmp\aidd-review-deps;src'
& D:\conda\python.exe -m pytest tests/test_macrocycle_full_audit.py tests/test_macrocycle_identity_audit.py tests/test_macrocycle_blocks.py tests/test_block_exploration.py tests/test_mol2_rdkit.py tests/test_mol2_import.py tests/test_chemical_companion.py -q
& D:\conda\python.exe scripts/check_english.py
git diff --check
```
