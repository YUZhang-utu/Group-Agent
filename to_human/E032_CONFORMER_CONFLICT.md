# E032: diagnose split_0119 conformer conflict

The 2026-09-14 failure is a repeated molecule-name/conformer-index pair with
different raw content hashes, not the random primary ID collision fixed by
`bab0643`. Do not rename, skip or overwrite records before comparing them.

Run in the existing workstation Python environment:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}" python scripts/diagnose_conformer_conflict.py \
  --registry /mnt/local/hand/yuzhang/aidd/library-precompute-20260911/registry.sqlite3 \
  --source /mnt/local/hand/yuzhang/aidd/mc_data/split_0119.mol2 \
  --molecule 'c--L-dA-Lnme-Wnme-dL-VNMe-c' \
  --conformer-index 1
```

If the source file is nested, use its actual full path. Return the JSON output.
The script opens SQLite read-only, streams the relevant sources, and compares
registered peers plus same-file peers (which may have disappeared from the
registry when the failed source transaction rolled back). Record indices are
zero-based. Differences are bounded to 100 lines per pair. A matching ordered
topology hash is not proof of chemical or stereochemical identity.

Default registration still rejects conflicting indices. The optional recovery
mode below retains differing records with audited internal indices.

Local verification: 15 focused tests passed, including ID collision retry,
whole-source rollback, exact duplicate behavior, and read-only diagnosis of
both committed and same-file peers. Real split_0119 records remain to be checked
on the workstation.

## Confirmed source collision and preservation mode (2026-09-15)

The workstation diagnostic found two records in split_0119.mol2 at zero-based
indices 131568 and 131625, both named
`c--L-dA-Lnme-Wnme-dL-VNMe-c_conf1`. Each has 116 atoms and 118 bonds;
ordered topology hashes match but coordinates and raw content hashes differ.
No committed registry peer was found, consistent with whole-source rollback.
This does not establish aligned geometry equivalence or stereochemical identity.

Resume with the additional `--preserve-index-conflicts` option:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
set -o pipefail
python scripts/precompute_library_batch.py \
  --source /mnt/local/hand/yuzhang/aidd/mc_data \
  --old-artifacts "$PWD/data/e019_artifacts/catalog.json" \
  --old-chemical /mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json \
  --output /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --workers 20 --preserve-index-conflicts --run \
  2>&1 | tee -a /mnt/local/hand/yuzhang/aidd/library-precompute-20260911/run.log
```

- No original source file, existing ID or completed artifact is changed.
- Distinct same-label records get separate conformer IDs and negative internal
  indices (-1, -2, ... within a molecule), avoiding collisions with later
  nonnegative source labels. The original name, content hash, file path and
  zero-based source position remain intact.
- `warnings_json` records original index, internal index and conflicting ID.
  Per-source progress and Ready logs expose `preserved_index_conflicts` counts.
- Exact content duplicates still use the existing duplicate rule. Ordered
  topology conflicts still stop and roll back the entire source transaction.
- Artifacts and chemical companions use source positions and conformer IDs,
  so internal relabeling does not select the wrong source record.
- Molecule grouping remains source-name-based. This mode preserves observations;
  it does not certify their chemical identity or remove the need for stereo QC.
- Completed sources are reused; the interrupted file restarts registration from
  its beginning. Repeating the resume command retains the same committed IDs.

Updated local suite: 136 passed, 1 RDKit-only module skipped. The real Linux
split_0119 build and full-library completion are pending workstation execution.
