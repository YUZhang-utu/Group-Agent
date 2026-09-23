# Molecule-budget library search and paged MOL2 handoff

This is a new explicit budget path. The existing threshold funnel is unchanged.
The current catalog contains 8,318,351 molecules / 25,813,808 conformers according
to the user's workstation report; it is not a 100-million-molecule library.

From Group-Agent in the `aidd-workstation` environment:

```bash
git pull --ff-only origin feature/structure-guided-chat
export AIDD_BUDGET_DIR=/mnt/local/hand/yuzhang/aidd/e059-budget-20260923
bash scripts/run_e059_budget.sh "$AIDD_BUDGET_DIR" 2>&1 | tee "$AIDD_BUDGET_DIR.log"
```

Use a persistent terminal session for the long run. Repeating the same command
resumes verified completed chunks. Changed code, inputs or settings require a
new output directory. Neither the library nor the original MOL2 files are edited.

Defaults: union retrieval budget 1,000,000 **unique molecules**, per-template
retrieval up to 1,000,000 molecules, RRF retrieval union cap, all stored conformers
of the selected molecules, all 11 templates, 24 workers. FAISS uses standardized
60-dimensional USRCAT, not the separate pharmacophore pair index. Search depth
and nprobe increase if necessary to obtain enough unique molecules, and actual
values are reported. Index checksum verification and loading are part of startup.
Current-target ANN/active recall is not established. No same-day completion claim.

The new path bypasses empirical coarse/contact/composite eligibility gates.
Contact scores still use the existing scoring definition. Protein collision
checking remains after seed generation; the existing seed enumeration cap remains.
Occupancy is annotation-only. Pose search is limited to stored conformers and
generated rigid seeds, not an exhaustive flexible docking calculation.

Final ranking: best actual pose per molecule/template; within-template molecule
ranks; a priority tier containing the union of each template's first 100,000;
RRF k=60 within each tier. The complete reserve tier is also retained. Correlated
templates and uncalibrated scores can bias ranking. A rank is a docking priority,
not evidence of activity. Fewer than 100,000 eligible molecules produces an explicit
shortfall, never duplicated molecules to fill the page.

Outputs:

- `retrieval.json`: actual library size, molecule count, search depths and probes.
- `status.json`: current completion/failure state.
- `chunks/`: sealed, resumable expensive computations.
- `ranking.sqlite`: complete molecule ranking, per-template ranks and scored poses.
- `page-000001-100000/molecules.csv`: rank, unique ID, original molecule and record
  names, representative conformer, file paths and chemical QC status.
- `molecules.jsonl`: all saved template representatives, score components,
  assignments, transformations, occupancy, source locations and chemistry audit.
- `original/MOL-*.mol2`: original source records, hash checked, names/H preserved.
- `posed/MOL-*.mol2`: the same records with the representative rigid transform
  applied to all atoms. These are not protonation-prepared docking inputs.

Representative selection for MOL2 uses best within-template rank, with template
ID ties; it does not compare uncalibrated raw scores between templates. Original
and posed files represent the same molecule and must not count as two molecules.
QC flags are retained for review. Source/companion topology or hash mismatch
stops export. New companion builds additionally store H counts and reconstruction
audits. Existing companion files are not silently rewritten; source-backed export
avoids their missing-H reconstruction problem.

Next nonoverlapping page, without rescoring:

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m aidd_agent.budget_export \
  --batch /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --run "$AIDD_BUDGET_DIR" --output "$AIDD_BUDGET_DIR/page-100001-200000" \
  --start-rank 100001 --count 100000
```

Do not merge pages from different ranking versions by rank number. Beyond the
stored pool, enlarge retrieval in a new run and deduplicate against delivered IDs.

The inexpensive existing pilot/reference distribution comparison is available:

```bash
python -m aidd_agent.pilot_review \
  --pilot /mnt/local/hand/yuzhang/aidd/pilot-256-20260923-135254 \
  --batch /mnt/local/hand/yuzhang/aidd/library-precompute-20260911 \
  --output /mnt/local/hand/yuzhang/aidd/e059-score-review --recover 0
```

It compares candidate molecule-best scores to crystal self-control scores. The
candidate population was already threshold-selected; pose optimization protocols
differ. This comparison does not establish activity discrimination or enrichment.
