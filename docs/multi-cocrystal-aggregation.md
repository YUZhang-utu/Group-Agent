# Multi-cocrystal 3D-search aggregation before docking

Run every co-crystal query independently through the tiered retrieval and
Gaussian pipeline. This stage combines completed detailed refinement results;
it never averages coordinates or replaces their source runs.

## 1. Write the query plan

Create `multi-cocrystal-plan.json`:

```json
{
  "format": "aidd-multi-cocrystal-query-plan",
  "version": 1,
  "library_id": "LIB-AFA68EE6888C",
  "queries": [
    {
      "query_id": "8BJU:QT9:A:401",
      "receptor_id": "8BJU-prepared-v1",
      "site_id": "WEE1-ATP-site",
      "result": "/path/query-8BJU/refine/merged-scores.npz",
      "result_manifest": "/path/query-8BJU/refine/merged-scores.manifest.json",
      "metadata": {"pdb_id": "8BJU", "ccd_id": "QT9"}
    },
    {
      "query_id": "OTHER:LIG:A:501",
      "receptor_id": "OTHER-prepared-v1",
      "site_id": "WEE1-ATP-site",
      "result": "/path/query-other/refine/merged-scores.npz",
      "result_manifest": "/path/query-other/refine/merged-scores.manifest.json",
      "metadata": {"pdb_id": "OTHER", "ccd_id": "LIG"}
    }
  ]
}
```

Use the same `site_id` only for biologically equivalent pockets. ATP and
allosteric pockets must have different site IDs. `receptor_id` identifies the
exact prepared receptor that will receive the docking task.

## 2. Aggregate and create docking admission

```bash
python -m aidd_agent.cli aggregate-multi-cocrystal-search \
  --plan /path/multi-cocrystal-plan.json \
  --output-dir /path/predocking/multi-cocrystal-v1 \
  --primary-objective atomcentered_anchored_joint \
  --top-conformers-per-query 3 \
  --per-query-quota 1000 \
  --consensus-quota 5000 \
  --global-limit-per-site 20000 \
  --rrf-k 60
```

The per-site global limit must be at least `number of queries in that site *
per-query-quota`. This makes query protection explicit instead of silently
dropping a co-crystal's unique chemical space.

## Outputs

- `query-molecule-evidence.jsonl`: conformer-to-molecule collapse for each
  query, independent best conformer/transform for every objective, and Top-M
  primary conformers;
- `molecule-summary.jsonl`: union across queries within each site, support
  counts, best rank, rank percentile and RRF score;
- `docking-admission.jsonl`: unique molecule admissions with protected-query
  and/or consensus reasons;
- `docking-tasks.jsonl`: one task per admitted molecule and supporting query,
  retaining exact receptor, site, conformer and candidate-to-query transform;
- `manifest.json`: hashes, parameters, counts and invariant checks.

Raw Gaussian scores remain available but are not averaged across co-crystal
queries. The cross-query order uses within-query molecule ranks, so chemically
different query ligands and score scales remain comparable. RRF is a reversible
priority view; all molecule evidence survives outside the admission budget.

