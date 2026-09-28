# GitHub checkout, full-library execution and mandatory validation

Repository: https://github.com/YUZhang-utu/Group-Agent

Branch: `feature/structure-guided-chat`.

Local release checks: 613 tests passed, 3 skipped; 56 focused macrocycle tests
passed. All three real mixed-source variants passed the new validation gates.
The small-source receipts are preserved in `E081_LOCAL_VALIDATION.json`.

Source root confirmed by the user:
`/mnt/local/hand/yuzhang/aidd/mc_data`. Available CSVs have the same stem and
directory as their MOL2. Missing CSVs retain explicit MOL2-only verification.

## Update an existing workstation checkout

Run from the AIDD repository root, preserving local work. If Git reports local
conflicts or divergence, stop and reconcile them; do not reset or force-pull.

```bash
git fetch origin
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
python -m pip install -e . --no-deps
python -m pytest tests/test_conformer_block_validate.py tests/test_verified_conformer_blocks.py tests/test_macrocycle_preflight.py tests/test_macrocycle_full_audit.py tests/test_macrocycle_identity_audit.py tests/test_macrocycle_blocks.py tests/test_block_exploration.py tests/test_mol2_rdkit.py tests/test_mol2_import.py tests/test_chemical_companion.py -q
```

Use the existing scientific environment with RDKit, NumPy and pytest. Record
`git rev-parse HEAD` with the run. `AIDD_PYTHON` can select a specific interpreter
for the shell runners; install and test with that same interpreter.

## One command for audit, blocks and validation

Choose a NEW run root outside the source tree. The example below is fixed to the
confirmed source path and does not run screening or docking:

```bash
bash scripts/run_e081_checked_library.sh /mnt/local/hand/yuzhang/aidd/mc_runs/e081-001
```

This runs source inventory and full audit, then backbone/chemistry/typed variants
at capacities 10000/20000/30000. Each variant must pass validation before the next
starts. No descriptors or memberships are silently overwritten. If interrupted,
use a new run root; resume scheduling is not implemented.

To use actual production IDs, append the registry SQLite path and library ID:

```bash
bash scripts/run_e081_checked_library.sh /mnt/local/hand/yuzhang/aidd/mc_runs/e081-002 /path/to/registry.sqlite3 ACTUAL_LIBRARY_ID
```

Without these arguments, IDs are explicitly source-only and cannot be treated as
production registry IDs. The supplied source root is known; the actual registry
path and library ID are not yet known. Registry matching is read-only and exact.

## Validation artifacts to inspect

Relative to the run root:

| Artifact | Meaning |
| --- | --- |
| `audit-stage/inventory/report.json` | File coverage, missing/unpaired CSVs and SHA-256 hashes |
| `audit-stage/audit/report.json` | All-source chemistry, stereo, name and conformer accounting |
| `blocks-stage/validation-typed/report.json` | Mandatory structural and release gates for all typed capacities |
| `blocks-stage/validation-typed/source_inspection.jsonl` | Names, source record indices, atom IDs, independently recomputed omega angles and flags for manual review |
| `blocks-stage/typed/review.jsonl` | Records withheld from fitted blocks and their reasons |
| `blocks-stage/typed/capacity-10000/memberships.jsonl` | Exact conformer/molecule IDs and block membership with provenance |

The same validation artifacts exist for backbone and chemistry variants. Require:

- `status: complete` and `structural_gate: passed`.
- `release_gate: passed_for_offline_use` for unattended continuation.
- The intended source manifest is complete; discovered files cannot prove that
  an omitted directory belongs to the intended library.
- Understand every missing conformer, boundary omega, symmetry and alignment flag.

`review_required` stops the wrapper with exit code 2 even if the structural checks
pass. Inspect and resolve the withheld records; do not reinterpret that receipt as
full acceptance. `hold` means validation failed. `errors` identifies the failure.
Do not hand-edit memberships or receipt hashes to force acceptance.

## What is checked to prevent wrong assignments

For EVERY extracted conformer and EVERY capacity:

1. Source admission identity, exact one-to-one accounting across descriptors and
   review, molecule/conformer IDs and source provenance.
2. Descriptor equality between extraction and fitted models; no extra or missing
   records and exact equality of exported memberships to the database.
3. Hard-stratum isolation, valid leaf assignment, capacity limits, complete tree
   populations and independent traversal of frozen split rules.
4. Recomputed leaf centroids and maximum distances, matching stored prototypes
   and radii; input/output hashes and receipt consistency.

For up to 200 deterministically selected source conformers per variant:

5. Reread original MOL2 records, check their identity/hash and atom maps, and
   reproduce the descriptor and residue/typed-feature mapping.
6. Independently compute omega using plane normals and compare angles and
   cis/trans/boundary states with the descriptor's projected-vector formula.

The report gives `source_samples_checked` and `raw_source_sampled_blocks` per
capacity. A 200-record inspection may not touch every block. It is not a claim
that every raw-source descriptor was independently recomputed. The full initial
audit and exhaustive membership checks have separate scopes. Inspect more records
when needed, with a fresh validation output:

```bash
python -m aidd_agent.conformer_block_validate \
  --build /mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/blocks-stage/typed \
  --output /mnt/local/hand/yuzhang/aidd/mc_runs/e081-typed-review-2000 \
  --source-samples 2000
```

## Integrity does not prove scientific clustering quality

Passing validates identity, bookkeeping, geometry reproduction and compliance with
the chosen partition rules. It does not establish optimal boundaries, energy
basins, biological activity or improved top-molecule recall. Different conformers
of one molecule are expected to occupy different blocks.

Before using blocks to skip expensive screening, follow the E080 replay protocol:
compare three variants/capacities on a held-out completed score reference, sample
100/500 unique molecules per block, use equal actual conformer-scoring budgets,
report both molecule discovery and ranking recall, and repeat seeds. No live
block rejection is enabled by this release. High maximum block radius and heavy
fragmentation are diagnostics for this comparison, not automatically calibrated
failure thresholds.

Validation adds full descriptor/database scans and sample source I/O. Plan disk
space for all three variants and three capacities. Full-workstation runtime and
memory have not been measured on this host.
