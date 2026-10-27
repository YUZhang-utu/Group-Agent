# KRAS: one million molecules per target for docking

The direct-count output is `/home/y/yuzhang/kras-direct-count-v1`. It stores
target/molecule IDs, not passing conformer IDs or poses. E118 must regenerate
bounded poses for conformers in that pool. This costs real compute; it is not
instant rescoring of cached poses and same-day completion is not guaranteed.

Original source results and reference packages are read-only. New scripts do not
modify the original E116 scorer or its input/code seals.

## Run the pilot first

From the repository with `aidd_workstation` activated:

```bash
python scripts/kras_soft_rank.py \
  --direct "$HOME/kras-direct-count-v1" \
  --output "$HOME/kras-soft-v1/pilot" \
  --sample 10000 --workers 20
cat "$HOME/kras-soft-v1/pilot/report.json"
```

Inspect each target's `eligible_molecules`, `geometry_poses`, score ranges and
the rough `full_hours_linear_estimate`. A zero-yield target blocks full execution.
This estimate cannot guarantee completion time or extrapolate distinct-molecule
yield. Native controls establish self recovery only, not activity enrichment.

## Full pool scoring, then export

```bash
nohup bash scripts/run_kras_million.sh \
  "$HOME/kras-direct-count-v1" "$HOME/kras-soft-v1" \
  > "$HOME/kras-soft-v1.log" 2>&1 < /dev/null &
```

The wrapper requires the completed pilot at `pilot/report.json`, resumes scoring
at `full/`, and exports to `export/`. Do not update code, original references or
libraries while a job runs. Scoring is resumable at atomic chunk boundaries.
Export uses a new directory and is not resumable: after an interrupted export,
preserve it and use the explicit export command with a different output path.

```bash
python scripts/export_kras_soft_rank.py \
  --ranking "$HOME/kras-soft-v1/full" \
  --output "$HOME/kras-soft-v1/export-retry" \
  --count 1000000 --batch-size 100000
```

Expected output when each target has at least one million eligible molecules:

```text
export/
  report.json
  selection.sqlite
  rejected.csv
  KRAS_G12D/
    batch_001.mol2 ... batch_010.mol2
    molecules.csv
  KRAS_G12V/
    batch_001.mol2 ... batch_010.mol2
    molecules.csv
```

Each target contains distinct library molecule IDs, one scored conformer per ID.
Shared molecules may occur in both target folders, with different representatives.
Files have 100,000 records except a final partial file if there is a shortfall.
`quota_met` and per-target `shortfall` disclose insufficient geometry survivors
or source-chemistry failures. A complete export does not necessarily meet quota.

## Scientific interpretation

- Original same-pose direct-contact gates remain, including mutant-specific rules.
- A 2.0 A heavy-atom minimum distance excludes gross overlaps; no shape/color hard
  cutoff and no hard reference-matching or regional-coverage cutoff remains.
- Same-pose score is an equal-block mean of shape, residue-group directional
  matching and mean reference-region coverage. It is an exploratory heuristic.
- Color is diagnostic only; reference-region coverage is not pocket volume.
- Retain three best conformers per molecule, exporting the highest scoring one.
- Select 80% by score, 20% by round-robin coarse graph similarity buckets. This
  diversity proxy is not a scaffold clustering or an activity ranking.
- Original MOL2 coordinates, atom/bond data and charges are preserved. Only the
  title changes. CSV records contain scores, source hashes, pose transforms and
  retained conformer IDs. Coordinates in MOL2 are NOT receptor-aligned; the saved
  transform records the evaluated pose. Apply docking-engine ligand preparation
  separately. No new protonation, minimization or docking has been performed.
- Raw-source access can be I/O intensive. Export reports source progress. Disk
  requirements depend on the selected molecular sizes and source filesystem.

## Checks and continuation

```bash
tail -n 30 "$HOME/kras-soft-v1.log"
cat "$HOME/kras-soft-v1/full/progress.json"
cat "$HOME/kras-soft-v1/export/report.json"
```

Return the pilot report before interpreting yield. At completion return both full
and export reports; keep molecule CSVs and MOL2 batches on the workstation for
docking. Do not transfer million-molecule files through GitHub.
