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

This change improves diagnosis only; it deliberately preserves the existing
conflict check. Completed batch sources are retained by the existing transaction
logic. Full batch resumption should wait for diagnosis and an appropriate fix.

Local verification: 15 focused tests passed, including ID collision retry,
whole-source rollback, exact duplicate behavior, and read-only diagnosis of
both committed and same-file peers. Real split_0119 records remain to be checked
on the workstation.
