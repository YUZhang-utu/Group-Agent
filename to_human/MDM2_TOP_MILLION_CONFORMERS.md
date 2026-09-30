# Select the top million saved MDM2 conformers

The user requests approximately one million CONFORMERS from the completed
3,854,146-conformer export. Do not select one million molecules, take the first
million source IDs, or use molecule ranking as a conformer score.

## Ranking evidence and policy

Historical ranking.sqlite compresses each molecule/template to one winning
conformer. Instead read search/chunks/TT/*.poses.jsonl before that compression.
The accompanying NPZ schedules and receipts cover every selected global ID for
every template. Validate file hashes, global/conformer/molecule identities,
stage-7 scored records, and full schedule coverage before ranking. If files are
missing or scored conformers are insufficient, stop without unscored padding.

Within each template rank unique conformers by optional_score descending,
gaussian_same_pose descending, then global ID ascending. Fuse these ranks with
RRF using the saved protocol's rrf_k and template_quota. Preserve the historical
positive-contact tier, then quota tier, then RRF ordering, now at conformer level.
Tie-break with global ID. This is a newly reconstructed conformer ranking, not
the old molecule ranking and not a claim of exhaustive library pose scoring.
Templates may have different scored populations. Missing scores contribute no
RRF term; count conformers without any score separately. Representative-template
contact/Gaussian columns describe the best within-template rank, not an affinity.

## Standalone execution

Download scripts/select_top_conformers.py and scripts/export_selected_conformers.py
to the same standalone directory; the latter supplies read/hash utilities only.
Dependencies: Python >=3.9, NumPy, SQLite with window functions and UPSERT support.
No RDKit, FAISS, original library rebuild, docking, or rescoring is performed.

```bash
SEARCH=/mnt/local/hand/yuzhang/aidd/e054-chat-20260923-104806/users/workstation/projects/prj-42c81b1322ab-prompt-aidd/runs/PROMPT-59addb9ffd6444c3/execution/guided/guided/search
EXPORTED=/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-1m-all-conformers
OUT=/mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers
python -u /mnt/local/hand/yuzhang/aidd/export-tools/select_top_conformers.py \
  --search "$SEARCH" --exported "$EXPORTED" --output "$OUT" --count 1000000
```

The output must be new. After interruption repeat with --resume and identical
code/inputs. Completed template aggregations and output parts are reused; an
interrupted template is reread. Do not run concurrent jobs against the same output.
The job reads many score chunks and exported MOL2 data; no runtime promise is
made. It needs temporary SQLite sort space and room for output MOL2. Source
directories and the running E094 checkout are not modified.

## Outputs and acceptance

- selection.sqlite: conformer-level aggregate scores, selected ranks and receipts
  for template processing. No original molecule ranking is overwritten.
- selected-conformers.csv: all selected IDs, rank, RRF, support, representative
  score evidence, original hash, output file and byte locator.
- top-00001.mol2 through top-00100.mol2 at default 10,000 conformers/file.
- Matching per-file CSV and hash receipt; original coordinates and names retained.
- report.json: input/scored/unscored/exported conformer counts, distinct molecule
  count, ranking policy and output hashes.

Verify status complete and exported_conformers 1000000; distinct_molecules may
be smaller. Output is original library geometry, not aligned or docked poses.
This ranking is restricted to the historical ANN-selected set, not the complete
8.3-million-molecule library. No current-target recall or affinity claim.

Local checks: seven tests passed across exporter and selector, including multiple
winning conformers from the same molecule despite opposite molecule ranks,
template-resume without double accumulation, missing score data, tampering, and
insufficient scored candidates. Full workstation execution remains pending.
