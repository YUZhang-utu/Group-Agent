# Review completed EquiScore leading blocks

Copy scripts/review_equiscore_leaders.py to the workstation repository scripts directory. It uses only the Python standard library and does not import or change scoring code, model environments or saved analyses.

Run from the workstation Group-Agent repository:

```bash
python3 scripts/review_equiscore_leaders.py --blocks 10 --output "$PWD/equiscore-leaders-review-20261006"
```

The default analysis is the existing equiscore-v3/equiscore-analysis directory from the handoff. Override with --analysis if the run was stored elsewhere. Use a different fresh output name if a review already exists; nothing is overwritten.

Outputs: review.md, leading_blocks.csv, leading_block_candidates.csv, report.json and a sibling ZIP bundle. Return the ZIP for interpretation. With nine eligible scheme/receptor groups and at least ten eligible blocks each, this selects 90 block records and 900 candidate rows; these are not necessarily 900 globally unique molecules. The primary block statistic remains the mean of the ten best distinct molecules. --blocks 5 changes only the number of displayed leading blocks. All 10 molecules remain exported, with ranks 1-5 flagged inspect_first.

The report lists source-sealed rankings, ChemPLP rank changes and candidate overlap. Original pose_file/pose_index values are retained without inferring the index base or exporting physical structures. Pose geometry, clashes, ligand chemistry and experimental performance require subsequent checks. Source reports or CSV hash mismatches stop export and must be investigated rather than bypassed. A partial full run is rejected. The analysis inputs are opened read-only.
