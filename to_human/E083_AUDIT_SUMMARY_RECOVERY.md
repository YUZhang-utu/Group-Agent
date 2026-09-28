# Recover a full audit interrupted during summary

The reported traceback is in the ring-size histogram, after ingestion commits.
It does not demonstrate incorrect conformer blocks: block construction has not
started in this invocation. SQLite GROUP BY may use temporary files on a different
filesystem from the database. Free space measured after a failed process exits
does not rule out transient exhaustion, a different temporary path or quotas.
Reference: https://www.sqlite.org/tempfiles.html

The user reports ample current free bytes and inodes on /tmp and the repository
NFS mount. The actual audit directory, environment and quota remain to be checked.
Do not delete the original audit or relaunch the full-library wrapper over it.
Do not run recovery while another process is writing to that audit or its sources.

## Get the repair and diagnose actual storage

From the repository root, with the existing AIDD Python environment active:

```bash
git pull --ff-only origin feature/structure-guided-chat
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"

# Replace this with the actual directory containing report.json and audit.sqlite.
AUDIT=/absolute/path/to/existing/audit
df -h "$AUDIT" /tmp
df -i "$AUDIT" /tmp
printenv SQLITE_TMPDIR TMPDIR TMP TEMP || true
quota -s
ulimit -f
ls -lh "$AUDIT/audit.sqlite" "$AUDIT/report.json" "$AUDIT/issues.jsonl"
```

NFS user/project quotas may require storage administrator confirmation even when
df shows available capacity. If environment variables select another filesystem,
check that path too. Do not infer that the source, repository and audit share a disk.

## Recover only the summary

Use a private local scratch directory on a filesystem with confirmed quota and
capacity. Set both variables before launching Python; this also applies to later
readiness/block stages that still use disk-backed temporary queries.

```bash
SQLITE_SCRATCH=$(mktemp -d /tmp/aidd-sqlite.XXXXXXXX)
export SQLITE_TMPDIR="$SQLITE_SCRATCH"
export TMPDIR="$SQLITE_SCRATCH"
python -m aidd_agent.macrocycle_full_audit --recover-summary --output "$AUDIT"
```

Recovery requires a saved failed report and complete source coverage. It checks
every source hash, SQLite quick_check, CSV occurrence totals, per-source conformer
counts and failure-report counters. It reads the existing database without
modification and does not repeat RDKit parsing/geometry work. Hashing still reads
all source bytes and the database, so recovery is not instantaneous.

Original report/issues are retained in a uniquely named before-summary-recovery
directory. Derived issue entries are rebuilt without duplicates. The complete
report is published only after successful summary and output hashing. Original
scientific code hashes remain intact; summary repair provenance is separate.
Partial ingestion, changed sources and completed reports are refused. If no failed
receipt could be written due to the storage failure, do not edit its status by hand.

The repair replaces three global JSON histogram sorts with streaming counters;
other exact duplicate/group checks can still need temporary disk. It does not
force all SQLite temporary data into RAM or claim to remove all disk requirements.

## Continue through the existing validation gates

If the original run used the E079/E081 wrapper, first compare inventory and audit:

```bash
python - "$AUDIT" <<'PY'
import json
import sys
from pathlib import Path
a = Path(sys.argv[1])
inventory = json.loads((a.parent / 'inventory/report.json').read_text())
audit = json.loads((a / 'report.json').read_text())
expected = {r['path']: r['sha256'] for r in inventory['files']}
observed = {r['path']: r['sha256'] for r in audit['csv_sources'] + audit['mol2_sources']}
assert audit['status'] == 'complete' and observed == expected, 'Inventory mismatch'
PY
```

Then use new, nonexistent output paths outside the source/audit trees:

```bash
python -m aidd_agent.macrocycle_preflight --audit "$AUDIT" \
  --output /absolute/new/readiness-recovered
bash scripts/run_e080_macrocycle_blocks.sh "$AUDIT" /absolute/new/blocks-recovered
```

Append REGISTRY_SQLITE and LIBRARY_ID to the final command if the original job
used a verified production registry. Each variant still runs mandatory E081 raw
source/assignment validation. A review_required result must be inspected; do not
bypass it. Recovery does not establish correct clustering or biological utility.
