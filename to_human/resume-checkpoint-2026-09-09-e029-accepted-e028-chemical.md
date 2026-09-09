# E029 accepted; E028 chemistry-aware QC continuation

Date: 2026-09-09

## Accepted E029 result

- 299,999 conformers across two source shards.
- Build time: 226.737 + 244.290 = 471.027 seconds (636.9 conformers/s).
- Exact source/output hashes, contiguous global IDs, boundary reads, and zero
  partials passed.
- Identical-command reuse: 1.49 seconds, 43,008 KB peak RSS, zero filesystem
  input, and unchanged catalog/shard manifest hashes.
- Catalog: `/mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json`.

The first-build peak RSS and aggregate CPU utilization were not captured, so
they are not claimed.

## Accepted E028 machine validation

- Output: `/mnt/local/hand/yuzhang/aidd/e028-predocking-qc-chemical-v1`.
- 200 selected tasks and 200 chemical SDFs: 100 each for 1X8B/824 and 8BJU/QT9.
- Chemistry/topology and vdW metrics are present; hashes passed.
- Retrieval membership and admission order remained unchanged.
- Close receptor points: 200/200; severe center-distance points: 195/200.
- Median candidate/query centroid displacement: 1.001 angstrom.

The collision counts do not invalidate retrieval. They establish that rigid
Gaussian overlays are hypotheses requiring local pose relaxation, not docking
inputs.

## Immediate continuation

Summarize vdW metrics by query and identify representative lowest, median, and
highest soft-exclusion poses:

```bash
python - <<'PY'
import json
from pathlib import Path
from statistics import median

root = Path('/mnt/local/hand/yuzhang/aidd/e028-predocking-qc-chemical-v1')
rows = [json.loads(line) for line in (root / 'pose-qc.jsonl').read_text().splitlines()]

for query in sorted({row['query_id'] for row in rows}):
    subset = [row for row in rows if row['query_id'] == query]
    penalties = [row['vdw_exclusion']['soft_exclusion_penalty'] for row in subset]
    fractions = [row['vdw_exclusion']['clashing_candidate_fraction'] for row in subset]
    severe = [row['vdw_exclusion']['severe_candidate_atoms'] for row in subset]
    ordered = sorted(subset, key=lambda row: row['vdw_exclusion']['soft_exclusion_penalty'])
    print('\nQUERY', query, 'n=', len(subset))
    print('penalty min/median/max=', min(penalties), median(penalties), max(penalties))
    print('clash_fraction median=', median(fractions),
          'severe_atoms median=', median(severe))
    for label, row in [('LOW', ordered[0]), ('MEDIAN', ordered[len(ordered)//2]),
                       ('HIGH', ordered[-1])]:
        print(label, 'rank=', row['query_molecule_rank'], 'gid=', row['global_id'],
              'penalty=', row['vdw_exclusion']['soft_exclusion_penalty'],
              'sdf=', row['chemical_sdf'])
PY
```

Inspect the LOW/MEDIAN/HIGH SDFs against their corresponding receptor in PyMOL.
Record whether collisions are local terminal-group conflicts, rigid-body frame
errors, or pervasive scaffold penetration. Only then integrate projected
direction and bounded terminal-torsion scoring into the fixed-budget final
refinement stage.

## Observed vdW distribution

- 1X8B/824: penalty min/median/max 13.942/50.197/83.300; median clashing atom
  fraction 0.582; median severe atoms 25.5.
- 8BJU/QT9: penalty min/median/max 6.900/17.337/49.025; median clashing atom
  fraction 0.389; median severe atoms 15.0.
- The high representatives have minimum center distances of 0.432 and 0.212
  angstrom, so at least those poses contain indisputable physical penetration.
- Six representative global IDs for visual review are 205388, 178819, 162174,
  269904, 297144, and 65851.
