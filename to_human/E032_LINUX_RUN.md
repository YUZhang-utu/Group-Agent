# E032 Linux run — expanded MOL2 library

Prepared 2026-09-11. **Not yet run on the real Linux library.**
User inventory: 365 MOL2, 277.08 GiB, 24 CPUs, approximately 58 GiB available RAM.

## Update the Linux checkout

The driver and controlled RDKit aromatic fallback are maintained in GitHub.
First preserve any manually copied untracked driver so it cannot block the pull:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git status --short
if [[ -f scripts/precompute_library_batch.py ]] && \
   ! git ls-files --error-unmatch scripts/precompute_library_batch.py >/dev/null 2>&1; then
  mv scripts/precompute_library_batch.py /tmp/precompute_library_batch.py.pre-e032b
fi
git pull --ff-only origin main
git rev-parse --short HEAD
```

Do not copy the old `/tmp` driver back. The fix includes the driver plus
`mol2.py`, `conformer_artifacts.py`, and `chemical_companion.py`. It first tries
full RDKit sanitization. Only if that fails does it run all operations except
Kekule assignment and preserve the original Tripos aromatic graph. Any other
sanitization failure still stops the shard. New manifests record
`rdkit_sanitization_counts`.

## Retry after the reported N5_0 failure

Repeat the pilot command below. Do not manually delete the partial shard or the
private registry. The driver moves `.N5_0.partial` into its owned `recovery/`
directory, reuses the transactional registry rows, and rebuilds that source.

## Pilot

Finish copying the MOL2 files before starting. Run on the compute workstation,
not an HPC login node. Use a persistent terminal such as an existing tmux session.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
conda activate aidd-workstation
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python scripts/precompute_library_batch.py \
  --source /mnt/local/hand/yuzhang/aidd/mc_data \
  --old-artifacts "$PWD/data/e019_artifacts/catalog.json" \
  --old-chemical /mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json \
  --output /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --workers 20 --max-new-files 1 --run
```

Do not pre-create the output directory: the driver establishes batch ownership.
The pilot verifies old artifacts and reconstructs the old stable identities by
reading the two original sources; it does not recompute their descriptors.
It then prepares one new file's v1 artifact, chemical companion and pharmacophore
postings. It does not train an index on an unrepresentative partial batch.

After it finishes, inspect both manifests:

```bash
python - <<'PY'
import json
from pathlib import Path
root = Path('/mnt/local/hand/yuzhang/aidd/library-precompute-20260911')
for path in (root/'artifacts/N5_0/manifest.json',
             root/'chemical/N5_0/manifest.json'):
    doc = json.loads(path.read_text())
    print(path, doc.get('conformers'), doc.get('rdkit_sanitization_counts'))
PY
```

Send the `Preserved incomplete output`, `Ready ...`, and `Pilot complete` lines,
plus both manifest-count lines. A completed pilot is not full-scale acceptance.

## Full run after a successful pilot

The same input inventory/output directory must be used. Omit `--max-new-files`.
This is the exact long-running command; `pipefail` preserves failure status.

```bash
set -o pipefail
python scripts/precompute_library_batch.py \
  --source /mnt/local/hand/yuzhang/aidd/mc_data \
  --old-artifacts "$PWD/data/e019_artifacts/catalog.json" \
  --old-chemical /mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json \
  --output /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --workers 20 --run \
  2>&1 | tee -a /mnt/local/hand/yuzhang/aidd/library-precompute-20260911/run.log
```

To resume after interruption, repeat the full command. Completed outputs are
hash-verified and reused. Incomplete directories created by this batch are moved
to its `recovery/` folder and rebuilt; no original data are deleted. Recovery is
per source/stage, not midway through a source. A POSIX lock rejects concurrent
writers. Do not add, replace or continue uploading inputs during the run.

## Outputs and boundaries

- `batch.json`: frozen input list and old catalog hashes.
- `registry.sqlite3`: private build registry; original deployment DB is untouched.
- `progress/*.json`: per-file record/duplicate counts and current-run elapsed time.
- `artifacts/catalog.json`: expanded catalog with unchanged old ID prefix.
- `chemical/catalog.json`: new and reused old E029 chemistry shards.
- `pharmacophore/catalog.json`: reusable target-independent pair postings.
- `faiss/trained/`: common transform and trained template sampled across all shards.
- `faiss/shards/`: independently encoded source shards.
- `faiss/index.faiss`, `transform.npz`, `manifest.json`: one combined final index.
- `COMPLETE.json`: end-to-end preprocessing completed; recall calibration pending.

Old v1/chemical shards remain at their original paths and are referenced, not
copied. Do not move or delete them. Pharmacophore postings are built in the new
snapshot; this driver does not assume the old pharmacophore catalog location.

The process can scan MOL2 more than once during initial registration, hashing,
artifact and companion construction. “One-time” means reusable preprocessing,
not literally one physical read pass. Finished query workflows use the binary
artifacts rather than parsing 277 GiB of text on every search.

Safety stops: changing inputs, repeated filename stems, ambiguous molecule /
conformer names, ordered topology conflicts, checksum or ID mismatch, fewer than
32 GiB free, or more than 1M new conformers in a single source (requires explicit
smaller shards). Nothing is silently renamed, chemically deduplicated or discarded.
Duplicates are only those recognised by the existing name/content identity rules.

The current registry is name-grouped; this pilot does not establish canonical
chemical identity across different names. Protonation, charge and stereochemical
quality remain scientific QC considerations for the imported library.

Build time is unknown until measured. Index building alone does not establish
recall or constant-time searches at larger scale. E031 remains query-dependent
and is not calculated for all molecules in advance. No docking/relaxation is run.

## Local verification

- Windows offline regression suite: 132 passed; one RDKit-only module was
  skipped in that interpreter.
- RDKit 2025.09.2 direct checks passed strict parsing, controlled aromatic
  fallback, and rejection of an independent invalid-valence record.
- A representative failing macrocycle passed both workers for all five
  conformers: 60 USRCAT values, 37 heavy atoms and 25 features per conformer.
- Linux multiprocessing, actual `N5_0.mol2`, FAISS training/merge and full-scale
  timing still require workstation validation.
