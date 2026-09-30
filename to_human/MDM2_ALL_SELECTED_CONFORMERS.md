# MDM2 selected molecule and conformer export

User-supplied workstation evidence on 2026-09-30 identifies two completed MDM2
budget searches, not exhaustive pose scoring. Both ranked 1,000,000 molecules
from the 8,318,351-molecule / 25,813,808-conformer indexed library and exported
the first 100,000 representatives. ANN recall on MDM2 remains unmeasured.

Workspace: `/mnt/local/hand/yuzhang/aidd/e054-chat-20260923-104806/users/workstation/projects/prj-42c81b1322ab-prompt-aidd/runs`.

- `PROMPT-59addb9ffd6444c3`: 11 templates; selection contains 3,854,146 conformers.
- `PROMPT-198eeebebde54547`: eight templates; all-conformer count not inspected.
- Batch: `/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`.
- Search relative path: `execution/guided/guided/search`.
- Existing representative export: `execution/guided/guided/page-000000001`.

`scripts/export_selected_conformers.py` is standalone, requiring Python >=3.9 and
NumPy. Copy/download it separately while E094 runs; do not update that running
checkout. Run once with `--batch BATCH --search SEARCH --output NEW_DIRECTORY`.
Use the same command with `--resume` after interruption. Never run concurrent
exporters into the same directory. A changed input or script prevents resume.

The exporter verifies saved selection/ranking hashes, shard ID file hashes,
global-ID to conformer/molecule identity, all registered conformer counts, and
every exported source record's ingestion-normalized SHA256. Missing original
source paths fail explicitly; files are not guessed by basename. No MOL2 parsing
or chemical reconstruction is needed. Original coordinates and names remain.
This does not establish chemical validity or docking pose quality.

Outputs: molecules.csv (all selected molecules, ranks, names, counts),
manifest.sqlite (full source mapping), part-NNNNN.mol2 (one per selected source
file), matching part CSVs (IDs, original source, record index, hash and output byte
offset/length), progress.json, per-part receipts and final report.json. Partitions
are storage files, not scientific conformer blocks. Text newline normalization is
the same as ingestion; output byte ranges recover hashed individual records.

Expected 11-template completion: status complete, molecules 1000000, conformers
3854146. The directory needs space for original MOL2 records and metadata.
This is additional sequential disk I/O; it can wait until E094 ends if resources
are shared. Completed parts are rehashed on resume, and incomplete parts rerun.

Local validation: three tests passed (all-conformer identity/count preservation,
completed-part tampering, missing conformer rejection, raw record tampering and
resume). English guard passed. Full workstation export has not run on this host.
