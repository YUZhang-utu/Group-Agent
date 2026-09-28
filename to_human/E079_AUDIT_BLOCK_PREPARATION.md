# Full-library audit and conformer block preparation

## Current evidence

E079 adds `aidd_agent.macrocycle_preflight`. It inventories source files and
prepares per-conformer admission from completed audit receipts. Production
clustering, registry joins and workstation full-library audit remain pending.
Earlier E077 evidence covers only the two supplied local shards, not the unseen
workstation library. No source data or production registry was modified.

## Workstation sequence

User confirmed on 2026-09-28: the complete workstation source directory is
`/mnt/local/hand/yuzhang/aidd/mc_data`; CSVs and MOL2s share the directory and stem.
This is user-provided coverage information, not a completed workstation scan.

After source delivery, run from the repository root:

```bash
bash scripts/run_e079_macrocycle_prepare.sh /mnt/local/hand/yuzhang/aidd/mc_audits/e079-run-001
```

The output root must not exist. Set `AIDD_PYTHON` to the scientific environment's
Python executable if needed. The script inventories sources, checks CSV headers,
runs the four-worker audit, compares every inventoried path/hash to audited inputs,
and generates readiness. It stops on failure and does not start clustering. The
script has not been executed on the workstation; local modules have fixture tests.

First deliver the reviewed working-tree source and tests. E075-E079 changes are
not committed or pushed; a workstation pull alone will not obtain them.
Use the repository scientific Python environment with this source installed.

1. Confirm the complete source root and place available CSVs adjacent to the
   matching MOL2, using the same stem and `.csv` suffix. Missing CSV is allowed.
2. Inventory files into a new directory outside the source tree:

```bash
python -m aidd_agent.macrocycle_preflight --source-dir /path/to/library --output /path/to/inventory-new
```

Review file counts, missing companions, unpaired CSVs and invalid CSV headers.
SHA-256 inventory reads every source byte but does not count conformers or perform
chemical validation. It does not prove that an omitted directory is covered.
Keep sources immutable during inventory and audit. Discovery currently expects
the companion extension `.csv`; differently cased extensions on case-sensitive
filesystems will be reported as unpaired, not silently matched.

3. Run the full audit into a separate new directory:

```bash
python -m aidd_agent.macrocycle_full_audit --source-dir /path/to/library --output /path/to/audit-new --workers 4
```

Compare every inventory source path/hash with `csv_sources` and `mol2_sources`
in the completed audit report. Resolve any added, removed or changed sources;
inventory is a receipt, not an enforced execution manifest. Review anomalies.
Interrupted audits still require fresh output; resumable scheduling is pending.

4. Generate block admission preparation:

```bash
python -m aidd_agent.macrocycle_preflight --audit /path/to/audit-new --output /path/to/readiness-new
```

This requires a completed audit and matching database/issues hashes. It streams
`conformer_readiness.jsonl`, preserving source path, record index, record hash,
original molecule/conformer names and peptide-ring source atom IDs. Duplicate
conformer names, same-name chemistry conflicts and failed audits remain review
records. Missing CSV alone does not reject geometry. Equivalent cyclic rotations,
boundary omega and unverified name alignment remain separate flags.

## Next implementation boundary

Join each readiness entry to production registry records using source provenance
and content hashes; do not use names alone. Translate source atom IDs into exact
artifact heavy-atom order and verify coordinates and chemical companion identity.
Then integrate the graph-verified peptide map into descriptor extraction. Retain
allowed cyclic rotations; graph-only names do not establish residue alignment.

Fit and benchmark actual conformer blocks at capacities 10000, 20000 and 30000.
The E079 report contains only optimistic block-count lower bounds ignoring strata,
not cluster memberships. Small residual blocks are allowed. Backbone, chemistry
and typed side-chain variants remain separate benchmarks. Molecule-aware block
sampling at 100/500 molecules must handle molecules appearing in multiple blocks;
evaluate molecule-level recall against completed screening, before block rejection.

## Local validation, 2026-09-28

34 focused tests passed, covering the audit, readiness, identity, block pilot,
exploration and MOL2/chemical adapters. Cross-record conflict tests prove that
individually accepted rows cannot bypass duplicate/chemistry review. Stale output
hashes and incomplete audit receipts are rejected before creating output.

Exploratory E078 receipt reuse produced 24 geometry candidates, preserving three
boundary-omega flags, including both CSV-backed and MOL2-only records. Inventory
covered two MOL2 files and one CSV, correctly retaining one missing companion.
Evidence: `data/e079-mixed-readiness/report.json` and
`data/e079-mixed-inventory/report.json`. No repeated two-shard audit was run.
