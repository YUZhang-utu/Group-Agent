# E029 workstation partial recovery and real-build continuation

Date: 2026-09-09

## Diagnosis

`chemical-companion-v1/.split_0001.partial` is an incomplete atomic-build
directory. It is not a resumable shard. Preserve it for diagnosis, move it out
of the output root, and rebuild `split_0001`. A completed `split_0001` directory,
if present, is reused automatically.

The original implementation could create an empty partial before rejecting a
missing relocated source. This is fixed locally: source existence, whole-file
SHA-256, and registry/artifact identity are now preflight checks. The recovery
commands below are also safe with the earlier implementation.

## 1. Inspect without deleting

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent

SOURCE_ROOT=/mnt/local/hand/yuzhang/aidd/mc_data
OUTPUT_ROOT=/mnt/local/hand/yuzhang/aidd/chemical-companion-v1
PARTIAL="$OUTPUT_ROOT/.split_0001.partial"

ls -lah "$SOURCE_ROOT/split_0001.mol2" "$SOURCE_ROOT/split_0002.mol2"
du -sh "$PARTIAL"
find "$PARTIAL" -maxdepth 1 -type f -printf '%f\t%s bytes\n' | sort
if [[ -f "$PARTIAL/manifest.json" ]]; then
  python -m json.tool "$PARTIAL/manifest.json"
fi
```

If `du` or `find` says the partial is absent, skip the move in step 3.

## 2. Verify the immutable sources

```bash
sha256sum "$SOURCE_ROOT/split_0001.mol2" "$SOURCE_ROOT/split_0002.mol2"
```

Required output:

```text
bc729d337ba7065ebeccdbd7aa05f2985d4ecc778dd8b7949e1b57ada49a3b2b  split_0001.mol2
6458081fa50eb2676d818497403084775fa13fbc87cbcfbe6c7318f40141d1e7  split_0002.mol2
```

Do not proceed if either digest differs.

## 3. Preserve the partial and rebuild

The move is recoverable and intentionally avoids `rm`:

```bash
RECOVERY=/mnt/local/hand/yuzhang/aidd/e029-recovery-20260909
mkdir -p "$RECOVERY"
if [[ -d "$PARTIAL" ]]; then
  mv "$PARTIAL" "$RECOVERY/.split_0001.partial"
fi

conda activate aidd-workstation
python -m pip install -e .
set -o pipefail

/usr/bin/time -v bash scripts/build_e029_chemical_companion.sh \
  /mnt/medchem_taltio/wrk/yu_agent/Group-Agent \
  --source split_0001="$SOURCE_ROOT/split_0001.mol2" \
  --source split_0002="$SOURCE_ROOT/split_0002.mol2" \
  2>&1 | tee /mnt/local/hand/yuzhang/aidd/e029-build-20260909.log
```

Because the calling shell enables `pipefail`, a build error remains visible
through the `tee` pipeline when the command is run exactly as shown.

## 4. Acceptance checks

```bash
test -f "$OUTPUT_ROOT/catalog.json"
find "$OUTPUT_ROOT" -maxdepth 1 -name '*.partial' -print
python -m json.tool "$OUTPUT_ROOT/catalog.json"

python - <<'PY'
import json
from pathlib import Path

root = Path('/mnt/local/hand/yuzhang/aidd/chemical-companion-v1')
catalog = json.loads((root / 'catalog.json').read_text())
assert len(catalog['shards']) == 2, catalog['shards']
assert sum(int(row['conformers']) for row in catalog['shards']) == 299999
for row in catalog['shards']:
    manifest = json.loads((root / row['name'] / 'manifest.json').read_text())
    assert manifest['source_sha256'] in {
        'bc729d337ba7065ebeccdbd7aa05f2985d4ecc778dd8b7949e1b57ada49a3b2b',
        '6458081fa50eb2676d818497403084775fa13fbc87cbcfbe6c7318f40141d1e7',
    }
    print(row['name'], manifest['conformers'], manifest['wall_seconds'],
          manifest['counts'])
print('total_conformers=299999')
PY
```

The `find` command should print nothing. Then rerun the identical timed build:
both completed shards should be reused, zero shard partials should appear, and
the catalog should remain valid. Preserve the first and reuse-run logs for the
confirmatory E029 record.

## Evidence to return

Return the two `sha256sum` lines, the pre-move partial inventory, the final JSON
printed by the build, GNU time metrics (elapsed time, CPU percent, maximum RSS),
the acceptance-check output, and the identical-command reuse-run output.
