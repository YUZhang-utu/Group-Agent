# Billion-scale tiered search

This execution path keeps broad and precise search products separately while
preventing total library size from controlling Gaussian cost. It is an
architecture-ready path; it is not yet a measured claim of one-billion-vector
latency or recall.

## Data plane

Each immutable library shard produces capped, already-ranked query results in a
directory:

```text
faiss-shard-000123/
  global_ids.npy   # little-endian int64, shape (local_K,)
  scores.npy       # float32/float64, higher is better, shape (local_K,)
```

Rows must be ordered by score descending and stable global ID ascending for
ties. The arrays are memory mapped. The full reusable library artifacts and
pharmacophore postings stay shard-local; only fixed ranked prefixes enter the
online merge.

Keep separate query result directories for:

- `baseline`: FAISS/L1; this prefix is mandatory and cannot be removed;
- `strict`: precise pharmacophore additions;
- `balanced`: intermediate additions;
- `loose`: broad exploration additions.

The source files and hashes remain in the schedule manifest. Budgeting controls
compute allocation, not a claim that omitted molecules cannot bind.

## 1. Merge ranked FAISS shards

Repeat `--shard-result` once per immutable shard:

```bash
python -m aidd_agent.cli merge-ranked-search-shards \
  --shard-result /mnt/local/hand/yuzhang/aidd/query-QT9/faiss-shard-000000 \
  --shard-result /mnt/local/hand/yuzhang/aidd/query-QT9/faiss-shard-000001 \
  --top-k 100000 \
  --output /mnt/local/hand/yuzhang/aidd/query-QT9/baseline-top100k.npz
```

This is an exact k-way merge over the available shard prefixes. A shard must
emit at least the local depth required by the frozen recall benchmark; global
exactness cannot recover a hit that a shard never emitted.

## 2. Build the fixed Gaussian schedule

Every input uses the same ranked schema. Omit a pharmacophore channel only when
it was deliberately not run.

```bash
python -m aidd_agent.cli build-tiered-candidate-schedule \
  --baseline /mnt/local/hand/yuzhang/aidd/query-QT9/baseline-top100k.npz \
  --strict /mnt/local/hand/yuzhang/aidd/query-QT9/strict-ranked \
  --balanced /mnt/local/hand/yuzhang/aidd/query-QT9/balanced-ranked \
  --loose /mnt/local/hand/yuzhang/aidd/query-QT9/loose-ranked \
  --baseline-budget 100000 \
  --strict-budget 100000 \
  --balanced-budget 50000 \
  --loose-budget 10000 \
  --output /mnt/local/hand/yuzhang/aidd/query-QT9/gaussian-schedule.npz
```

The maximum Gaussian admission is 260,000 conformers with these defaults,
independent of whether the indexed library contains 10 million or one billion
conformers. Duplicates consume no addition budget and accumulate source flags.

## 3. Retain slim coarse search and detailed refinement

```bash
python -m aidd_agent.cli run-scaled-gaussian-reranking \
  --artifact-catalog /mnt/medchem_taltio/wrk/yu_agent/Group-Agent/data/e019_artifacts/catalog.json \
  --query /mnt/medchem_taltio/wrk/yu_agent/Group-Agent/data/e019_query_8bju/gaussian-query-v1.npz \
  --candidate-schedule /mnt/local/hand/yuzhang/aidd/query-QT9/gaussian-schedule.npz \
  --output-dir /mnt/local/hand/yuzhang/aidd/query-QT9/gaussian-scaled-v1 \
  --workers 16 \
  --coarse-chunk-size 2000 \
  --refine-chunk-size 250 \
  --top-n-per-objective 5000 \
  --max-pair-seeds 512
```

Coarse chunk rows contain one int64 ID and three float32 objective scores, a
20-byte logical payload per conformer. They remain on disk after refinement.
The streaming selector is exact relative to a full stable sort but does not
build a monolithic detailed coarse archive. Detailed transforms, seed counts,
and raw overlaps are retained in the refined result only.

Use `--stage coarse` and `--stage refine` for separate scheduler jobs. Repeating
the identical command reuses hash-valid chunks. Use a new output directory when
inputs or scientific parameters change.

## Production gates

Before calling the path billion-scale validated, run the frozen ladder:

1. reproduce E023 IDs and Top-N selection on 299,999 QT9 conformers;
2. measure merge time, peak RSS, bytes/row and resume at 10M and 100M;
3. validate per-shard FAISS local depth against exhaustive truth;
4. run at least one true multi-shard production query;
5. physically benchmark one billion vectors before reporting billion-scale
   latency or recall.

